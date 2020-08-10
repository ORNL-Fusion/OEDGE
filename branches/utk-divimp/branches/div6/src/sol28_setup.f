c     -*Fortran*-
c
c ======================================================================
c
      SUBROUTINE SetTargetConditions(itube)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      INTEGER,INTENT(IN) :: itube

      INTEGER i
      LOGICAL first_call
      DATA    first_call /.TRUE./

      tube(itube)%jsat(2,1:S28_MAXNION) = 0.0

      IF (first_call) THEN
        first_call = .FALSE.
c...    This is nasty, but it's necessary to initialize all of the targets,
c       otherwise any regions linked to other regions, where it's assumed the 
c       linked-to region has already been defined, will fail.  This is a backwards
c       compatibilty issue at the moment and only affects the legacy LPDATI,O 
c       assignment for SOL28, i.e. if doesn't bother DIVIMP or the newer SOL28
c       target data assignments:
        DO i = grid%isep, ntube
          CALL SetTargetConditions_Legacy(i)        
        ENDDO
      ENDIF

      CALL SetTargetConditions_Legacy(itube)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE AssignReferenceSolution
      USE mod_sol28
      USE mod_sol28_global
      IMPLICIT none

      INTEGER ion

c...  Copy data to REF_xxx variables:
      ref_ntube  = ntube
      ref_nion   = nion
      ref_ncell  = ncell  
      ref_nfluid = nfluid
      IF (ncell.NE.nfluid) 
     .  CALL ER('AssignReferenceSolution','NCELL.NE.NFLUID',*99)
      IF (ALLOCATED(ref_tube )) DEALLOCATE(ref_tube )
      IF (ALLOCATED(ref_cell )) DEALLOCATE(ref_cell ) 
      IF (ALLOCATED(ref_fluid)) DEALLOCATE(ref_fluid)
      ALLOCATE(ref_tube (ntube))
      ALLOCATE(ref_cell (ncell))
      ALLOCATE(ref_fluid(nfluid,nion))
      ref_tube(1:ntube) = tube(1:ntube)
      DO ion = 1, nion  
        ref_cell (1:ncell)      = cell (1:ncell)
        ref_fluid(1:nfluid,ion) = fluid(1:nfluid,ion)
      ENDDO

c      WRITE(0,*) 'REF:',fluid(1:100,1)%ne
c      WRITE(0,*) 'REF:',ref_fluid(1:100,1)%ne

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadReferenceSolution(mode)
      USE mod_sol28_params
      USE mod_sol28
      USE mod_sol28_global
      IMPLICIT none

      INTEGER, INTENT(IN) :: mode

      INTEGER        ion,status
      CHARACTER*1024 fname,command

      IF (opt%osm_load.EQ.0) THEN
        ref_ntube  = 1
        ref_nion   = 1
        ref_nfluid = 1
        ALLOCATE(ref_tube (ref_ntube))
        ALLOCATE(ref_fluid(ref_nfluid,ref_nion))
        RETURN
      ENDIF

c...  Store current OSM plasma:
      CALL SaveGrid('osm.tmp')

c...  Load reference OSM plasma (binary form):
      fname = TRIM(opt%f_osm_load)

c...  Copy file to run directory:
      command = 'cp '//TRIM(opt%f_osm_dir)//TRIM(fname)//' .'
      CALL CIssue(TRIM(command),status)
      IF (status.NE.0) 
     .  CALL ER('LoadReferencePlasma','Unable to copy plasma file',*99)

c...  Unzip file if necessary:
      IF (fname(LEN_TRIM(fname)-2:LEN_TRIM(fname)).EQ.'zip'.OR.
     .    fname(LEN_TRIM(fname)-1:LEN_TRIM(fname)).EQ.'gz') 
     .  CALL UnzipFile(fname)

      SELECTCASE (mode)
        CASE (1) 
c...      Load reference solution:
          WRITE(0,*) 'LOADING REFERENCE PLASMA:'      
          WRITE(0,*) '  FILE NAME = ',TRIM(fname)
          CALL LoadGrid(TRIM(fname))
          CALL AssignReferenceSolution
        CASE DEFAULT
          CALL ER('LoadReferenceSolution','Unknown MODE',*99)
      ENDSELECT

c...  Need to add a better clean-up of the reference solution somewhere,
c     since it probably isn't needed all the time and might take up
c     a lot of space... currently it is deallocated in CLEANUP...

c...  Restore current OSM grid:
      CALL LoadGrid('osm.tmp')

      RETURN
 99   WRITE(0,*) '  FILE NAME = ',TRIM(fname)
      WRITE(0,*) '  COMMAND   = ',TRIM(command)
      WRITE(0,*) '  ERROR     = ',status
      STOP
      END
c
c ======================================================================
c (from InterpolateTargetdata in div6/src/setup.f -SL, 09.01.07)
c
c
      SUBROUTINE SetTargetConditions_Legacy(itube)
      USE mod_legacy
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      INTEGER,INTENT(IN) :: itube

c.... OSM:
      INTEGER ishift


c...  Local:
      REAL GetRelaxationFraction

      INTEGER iitersol

      INTEGER ii,ind1,ind2,i1,i2,i3,region,method,ring,target,idat,
     .        no_data_warning,ir
      LOGICAL specific,apply,cycle_loop,initialise
      REAL    dum1,dum2,dum3,dum4,exp2,exp3,exp4,dpsin,psin0,frac,
     .        tedat(MAXNRS),tidat(MAXNRS),nedat(MAXNRS)

      DATA no_data_warning, initialise / 0, .TRUE. / 
      SAVE


      IF (itube.LT.grid%isep) RETURN


      CALL AssignLegacyVariables

      IF (logop.GT.0) WRITE(logfp,*) 'SETTING TARGET CONDITIONS'
      
      IF (tarninter(LO).EQ.0.OR.tarninter(HI).EQ.0) THEN
c        CALL ER('InterpolateTargetData','No data to '//
c     .          'interpolate?',*99)
c        CALL WN('InterpolateTargetData','No data to '//
c     .          'interpolate?')
        WRITE(logfp,*) 'WARNING InterpolateTargetData: No data '//
     .                 'to interpolate'
        RETURN
      ENDIF

      IF (connected) THEN
        WRITE(0,*) 'CHEATING ON PSITARG'
c        DO ir = irsep, nrs
          psitarg(ir,1) = REAL(ir-irsep)/100.0 + 1.0
          psitarg(ir,2) = REAL(ir-irsep)/100.0 + 1.0
c        ENDDO
      ENDIF

      WRITE(SLOUT,*)
      WRITE(SLOUT,*) 'INTERPOLATING TARGET DATA'

c...  Need to check that PSITARG is assigned:
      IF (psitarg(irsep,2).EQ.0.0)
     .  CALL WN('InterpolateTargetData','PSITARG does not appear to '//
     .                                  'be assigned')

c      repeat = .TRUE.

      frac = GetRelaxationFraction()

      IF (initialise) THEN
        WRITE(logfp,*) 'INITIALISING SetTargetConditions_Legacy'
        WRITE(logfp,*) '  IRSEP= ',irsep
        WRITE(logfp,*) '  NRS  = ',nrs
        ii = 0
        DO ir = irsep, nrs
          ii = ii + 1
          lpdati(ii,1) = REAL(ir)
          lpdati(ii,2) = 1.0
          lpdati(ii,3) = 1.0
          lpdato(ii,1) = REAL(ir)
          lpdato(ii,2) = 1.0
          lpdato(ii,3) = 1.0
          IF (lpdatsw.EQ.0) THEN
            lpdati(ii,4) = 1.0E+12
            lpdato(ii,4) = 1.0E+12
          ELSE
            lpdati(ii,4) = 1.0E+00
            lpdato(ii,4) = 1.0E+00
          ENDIF
        ENDDO
        nlpdato = ii
        nlpdati = ii
        initialise = .FALSE.
      ENDIF


c...  Assign target data (LPDATO, LPDATI) from the TARINTER
c     data listed in the DIVIMP/OEDGE input file:
      ir = tube(itube)%ir

 10   DO region = IKLO, IKHI
        IF (tarninter(region).EQ.0) CYCLE
        DO ii = 1, nlpdato
          IF (NINT(lpdato(ii,1)).EQ.tube(itube)%ir) EXIT
        ENDDO
        IF (ii.EQ.nlpdato+1) CALL ER('SetTargetConditions_Legacy',
     .                               'Ring not identified',*99)

c...    Check if target data has been assigned specifically for this ring: 
        specific = .FALSE.
        method = 1
        ind1 = 0
        ind2 = 0

        DO i1 = 1, tarninter(region)
          apply = .FALSE.
          IF (tarinter(i1,1,region).EQ.-1.0) THEN
            specific = .TRUE.

            IF     (ir.GE.NINT(tarinter(i1,2,region)).AND.
     .              ir.LE.NINT(tarinter(i1,3,region))) THEN
c...          Apply target data to ring:
              apply = .TRUE.
            ELSEIF (tarinter(i1,2,region).LT.0.0.AND.
     .              tarinter(i1,3,region).LT.0.0) THEN
c...          Make sure that ring group exists:
              IF (-NINT(tarinter(i1,3,region)).GT.grdntreg(region)) THEN
                CALL WN('InterpolateTargetData','Specified '//
     .                  'region index is not valid')
                RETURN
              ENDIF
c     .          CALL ER('InterpolateTargetData','Specified '//
c     .                  'region index is not valid',*99)
c...          Check if ring IR is a member of the specified group(s)
c             of rings:
              DO i2 = -NINT(tarinter(i1,2,region)), 
     .                -NINT(tarinter(i1,3,region))
                DO i3 = 1, grdntseg(i2,region)
                  IF (ir.EQ.grdtseg(i3,i2,region)) apply = .TRUE. 
                ENDDO
              ENDDO                           
            ENDIF
            IF (apply) THEN
c...          Interpolation scheme:
              IF (tarinter(i1,4,region).EQ.1.0) THEN              
c...            Versus PSIn, exponential decay:
                method = 2
              ELSE
c...            Standard:
                method = 1
              ENDIF
              ind1 = 0
              ind2 = tarninter(region)
              DO i2 = i1+1, tarninter(region)-1
                IF (ind1.EQ.0.AND.
     .              tarinter(i2  ,1,region).NE.-1.0) ind1 = i2
                IF (ind1.NE.0.AND.
     .              tarinter(i2+1,1,region).EQ.-1.0) THEN
                  ind2 = i2 
                  EXIT
                ENDIF
              ENDDO
              WRITE(PINOUT,*) 'APPLYING:',ir,ind1,ind2,region
              IF (ind2-ind1+1.LT.2) 
     .          CALL ER('InterpolateTargetData','Insufficient '//
     .                  'target data for interpolation',*99)
            ENDIF
          ELSEIF (tarinter(i1,1,region).EQ.-2.0.AND.
     .            ((NINT(tarinter(i1,5,region)).EQ.ir.AND.
     .              NINT(tarinter(i1,6,region)).EQ.0).OR.
     .             (NINT(tarinter(i1,5,region)).LE.ir.AND.
     .              NINT(tarinter(i1,6,region)).GE.ir))) THEN
c...        Assign data directly to the ring, without interpolation:
            method = 3
            apply = .TRUE.
            idat = i1
          ENDIF
        ENDDO

        IF (.NOT.specific) THEN  
c...      Target data was not specified for specific rings, so treat
c         the entire TARINTER array as target data to be interpolated:
          ind1 = 1
          ind2 = tarninter(region)
        ELSEIF (ind1.EQ.0.OR.ind2.EQ.0) THEN
c...      Data was not found for this ring:
          IF (no_data_warning.EQ.0) no_data_warning = 1
          WRITE(88,*) 'WARNING: TARGET DATA NOT FOUND FOR',ir,region
          CYCLE
        ENDIF

        IF     (method.EQ.1) THEN
c...      Linearly interpolate target data from TARINTER arrays:
          IF (region.EQ.IKLO) THEN
            dum1 = psitarg(ir,2)
           ELSE
            dum1 = psitarg(ir,1)
          ENDIF

c...      Make a list:            
          i3 = 0
          DO i2 = ind1, ind2
            i3 = i3 + 1
            tedat(i3) = (1.0 - frac) * tarinter(i2,2,region) + 
     .                         frac  * tarinter(i2,5,region)
            tidat(i3) = (1.0 - frac) * tarinter(i2,3,region) + 
     .                         frac  * tarinter(i2,6,region)
            nedat(i3) = (1.0 - frac) * tarinter(i2,4,region) + 
     .                         frac  * tarinter(i2,7,region)
          ENDDO
          CALL Fitter(ind2-ind1+1,tarinter(ind1,1,region),tedat,
     .                1,dum1,dum2,'LINEAR')
          CALL Fitter(ind2-ind1+1,tarinter(ind1,1,region),tidat,
     .                1,dum1,dum3,'LINEAR')
          CALL Fitter(ind2-ind1+1,tarinter(ind1,1,region),nedat,
     .                1,dum1,dum4,'LINEAR')
        ELSEIF (method.EQ.2) THEN
c...      Exponential decay along the target versus PSIn:
          target = 0
          IF (tarinter(ind1,1,region).LT.0.0) THEN
c...        Decay is scaled by the values for a specified target segement:
            ring   = -NINT(tarinter(ind1,1,region))
            target =  NINT(tarinter(ind1,2,region))
            IF     (target.EQ.1) THEN
              DO i2 = 1, nlpdato
                IF (NINT(lpdato(i2,1)).EQ.ring) THEN
                  tarinter(ind1,2,region) = lpdato(i2,2) / te_mult_o
                  tarinter(ind1,3,region) = lpdato(i2,3) / ti_mult_o
                  tarinter(ind1,4,region) = lpdato(i2,4) /  n_mult_o
                  tarinter(ind1,5,region) = tarinter(ind1,2,region)
                  tarinter(ind1,6,region) = tarinter(ind1,3,region)
                  tarinter(ind1,7,region) = tarinter(ind1,4,region)
                ENDIF
              ENDDO
              psin0 = psitarg(ring,2)
            ELSEIF (target.EQ.2) THEN
              DO i2 = 1, nlpdati
                IF (NINT(lpdati(i2,1)).EQ.ring) THEN
                  tarinter(ind1,2,region) = lpdati(i2,2) / te_mult_i
                  tarinter(ind1,3,region) = lpdati(i2,3) / ti_mult_i
                  tarinter(ind1,4,region) = lpdati(i2,4) /  n_mult_i
                  tarinter(ind1,5,region) = tarinter(ind1,2,region)
                  tarinter(ind1,6,region) = tarinter(ind1,3,region)
                  tarinter(ind1,7,region) = tarinter(ind1,4,region)
                ENDIF
              ENDDO
              psin0 = psitarg(ring,1)
            ELSE
              CALL ER('InterpolateTargetData','Invalid target',*99)
            ENDIF
          ELSE
            psin0 = tarinter(ind1,1,region)
          ENDIF
          IF (region.EQ.IKLO) THEN
            dpsin = ABS(psin0 - psitarg(ir,2))
          ELSE
            dpsin = ABS(psin0 - psitarg(ir,1))
          ENDIF
c...      Interpolate target profile parameters (in case target boundary
c         relaxation is being used):
          DO i2 = 1, 3
            tedat(i2) = (1.0 - frac) * tarinter(ind1+i2-1,2,region) +
     .                         frac  * tarinter(ind1+i2-1,5,region)
            tidat(i2) = (1.0 - frac) * tarinter(ind1+i2-1,3,region) +
     .                         frac  * tarinter(ind1+i2-1,6,region)
            nedat(i2) = (1.0 - frac) * tarinter(ind1+i2-1,4,region) +
     .                         frac  * tarinter(ind1+i2-1,7,region)
          ENDDO
c...      Allow over-ride (RELMODE=8) of profile parameters at each step,
c         when boundary condition relaxation is being used, REL_OPT=2,3: 
          IF (tedat(1).EQ.-1.0) THEN
            IF (relmode.EQ.8.AND.relexpdat(region,1,ir).NE.0.0) THEN
              tedat(1) = relexpdat(region,1,ir)
            ELSEIF (rel_step.EQ.0) THEN
c...          The 0th EIRENE iteration is (currently) being avoided when target relaxation
c             is being used, so just assign some arbitrary nedat(1)!=-1.0.  The value from 
c             RELEXPDAT will be assigned on a subsequent call to this routine before the 
c             plasma is calculated (when REL_STEP=1, see in BgPlasma routine):
              tedat(1) = 1.0E+01
            ELSE
              CALL ER('InterpolateTargetData','Bad Te over-ride',*99)
            ENDIF
          ENDIF
          IF (tidat(1).EQ.-1.0) THEN
            IF (relmode.EQ.8.AND.relexpdat(region,2,ir).NE.0.0) THEN
              tidat(1) = relexpdat(region,2,ir)
            ELSEIF (rel_step.EQ.0) THEN
              tidat(1) = 1.0E+01
            ELSE
              CALL ER('InterpolateTargetData','Bad Ti over-ride',*99)
            ENDIF
          ENDIF
          IF (nedat(1).EQ.-1.0) THEN
            IF (relmode.EQ.8.AND.relexpdat(region,3,ir).NE.0.0) THEN
              nedat(1) = relexpdat(region,3,ir)
            ELSEIF (rel_step.EQ.0) THEN
              nedat(1) = 1.0E+01
            ELSE
              CALL ER('InterpolateTargetData','Bad ne over-ride',*99)
            ENDIF
          ENDIF

c...      Set exponents:
          exp2 = EXP(-dpsin / tedat(3))
          exp3 = EXP(-dpsin / tidat(3))
          exp4 = EXP(-dpsin / nedat(3))
c...      Calculate target boundary conditions:
          dum2 = (tedat(1) - tedat(2)) * exp2 + tedat(2)
          dum3 = (tidat(1) - tidat(2)) * exp3 + tidat(2)
          dum4 = (nedat(1) - nedat(2)) * exp4 + nedat(2)
c...      (Need to replace the target data specifier which was overwritten 
c          above, nasty):
          IF (target.NE.0) tarinter(ind1,2,region) = REAL(target)
        ELSEIF (method.EQ.3) THEN
c...      Direct application of target values:
          dum2 = tarinter(idat,2,region) 
          dum3 = tarinter(idat,3,region) 
          dum4 = tarinter(idat,4,region) 
        ENDIF
c...    Assign LPDATx arrays (target data used to assign target data
c       KxDS arrays):
        IF (region.EQ.IKLO) THEN
          nlpdato = MAX(nlpdato,ii)
          lpdato(ii,1) = REAL(ir)
          lpdato(ii,2) = dum2 * te_mult_o
          lpdato(ii,3) = dum3 * ti_mult_o
          lpdato(ii,4) = dum4 *  n_mult_o
        ELSE
          nlpdati = MAX(nlpdati,ii)
          lpdati(ii,1) = REAL(ir)
          lpdati(ii,2) = dum2 * te_mult_i
          lpdati(ii,3) = dum3 * ti_mult_i
          lpdati(ii,4) = dum4 *  n_mult_i
        ENDIF
c...    End of REGION loop:
      ENDDO

c      IF (repeat) THEN
c...    Need to do this twice:
c        repeat = .FALSE.
c        GOTO 10
c      ENDIF

c      CALL Outputdata(85,'sdfds')
c      STOP 'sdgsdg'
	
c...  Copy data to OSM arrays:
c      DO itube = 1, ntube
c        IF (itube.LT.grid%isep) CYCLE

        cycle_loop = .TRUE.
        DO i1 = 1, opt%sol_n
          IF (itube.GE.opt%sol_tube(1,i1).AND.
     .        itube.LE.opt%sol_tube(2,i1)) cycle_loop = .FALSE.
        ENDDO
c        IF (cycle_loop) CYCLE

c        WRITE(0,*) 'WHAT?',itube,opt%sol_n
c        WRITE(0,*) 'WHAT?',i1,opt%sol_tube(1:2,1)
c        WRITE(0,*) 'WHAT?',i1,opt_iteration(1)%sol_tube(1:2,1)

        IF (.NOT.cycle_loop) THEN

c          DO ii = 1, nlpdato
c            ishift = 1
c            IF (itube.GE.grid%ipfz) ishift = 2
c             WRITE(0,*) 'LPDATI:',itube,tube(itube)%ir
c            IF (NINT(lpdato(ii,1)).EQ.tube(itube)%ir) THEN
c              WRITE(0,*) '      :',ii,lpdato(ii,2)
              tube(itube)%te  (LO)   = lpdato(ii,2)
              tube(itube)%ti  (LO,1) = lpdato(ii,3)
              tube(itube)%jsat(LO,1) = lpdato(ii,4)
              WRITE(logfp,*) 'JSAT8:',itube,tube(itube)%ir,lpdato(ii,4)
c            ENDIF
c          ENDDO

c          DO ii = 1, nlpdati
c            ishift = 1
c            IF (itube.GE.grid%ipfz) ishift = 2
c            IF (NINT(lpdati(ii,1)).EQ.tube(itube)%ir) THEN
c            IF (NINT(lpdati(ii,1)).EQ.itube+ishift) THEN
              tube(itube)%te  (HI)   = lpdati(ii,2)
              tube(itube)%ti  (HI,1) = lpdati(ii,3)
              tube(itube)%jsat(HI,1) = lpdati(ii,4)
              WRITE(logfp,*) 'JSAT7:',itube,tube(itube)%ir,lpdati(ii,4)
c            ENDIF
c          ENDDO

       ENDIF

c      ENDDO

c      STOP 'fgdfsg'

      IF (no_data_warning.EQ.1) THEN
        WRITE(0,*)
        WRITE(0,*) '****************************************'
        WRITE(0,*) '* TARGET DATA NOT FOUND FOR SOME TUBES *'
        WRITE(0,*) '****************************************'
        WRITE(0,*)
        no_data_warning = 2
      ENDIF

      RETURN
 99   WRITE(0,*) '  I1,REGION   =',i1,region
      WRITE(0,*) '  ITUBE       =',itube
      WRITE(0,*) '  IR          =',ir
      DO i1 = 1, nlpdato
        WRITE(0,*) 'NLPDATO:',i1,lpdato(i1,1)
      ENDDO
      DO i1 = 1, nlpdati
        WRITE(0,*) 'NLPDATO:',i1,lpdati(i1,1)
      ENDDO
 
c      WRITE(0,*) ' -NINT(TARINTER)  =',-NINT(tarinter(i1,3,region))
c      WRITE(0,*) '  GRDNTREG(REGION)=',grdntreg(region)
      STOP
      END
c
c ======================================================================
c 
c subroutine: FindCell_New
c
c Identifies the first, last and current cell that the node line segment
c intersects (I think) -- sampling from the stated range of tubes for
c the interpolation node, and not outside this range.
c
      SUBROUTINE FindCell_New(ind0,ind1,itgive,iccell)
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_geometry
      IMPLICIT none

      INTEGER ind0,ind1,itgive,iccell(3)

      INTEGER ic,it,i1,ic1(ntube),iobj,isrf,ivtx(2),itcell(3)
      REAL    dist,dist1(ntube),clcell(3)
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd

      ic1 = 0
      dist1 = 0.0

      DO it =     osmnode(ind1)%tube_range(1), 
     .        MIN(osmnode(ind1)%tube_range(2),ntube)
        IF (tube(it)%type.EQ.GRD_BOUNDARY) CYCLE
c        WRITE(0,*) 'FC: IT=',it
        DO i1 = ind0+1, ind1
c          WRITE(0,*) 'FC: I1=',i1
          a1 = DBLE(osmnode(i1-1)%rad_x)  ! Collect line segment end points from the node setup data
          a2 = DBLE(osmnode(i1-1)%rad_y)
          b1 = DBLE(osmnode(i1  )%rad_x)
          b2 = DBLE(osmnode(i1  )%rad_y)
          dist = REAL(DSQRT((a1 - b1)**2 + (a2 - b2)**2))
          DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)
c...        Assumed 1:1 mapping between grid and data:
            iobj = ic
            isrf = ABS(obj(iobj)%iside(1))
            ivtx(1:2) = srf(isrf)%ivtx(1:2)
            c1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))  ! *** Clean up with GetVertex calls ***
            c2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2)))
            isrf = ABS(obj(iobj)%iside(3))
            ivtx(1:2) = srf(isrf)%ivtx(1:2)
            d1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
            d2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2)))
            CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
            IF (tab.GE.0.0.AND.tab.LT.1.0.AND.
     .          tcd.GE.0.0.AND.tcd.LT.1.0) THEN
              ic1(it) = ic
              dist1(it) = dist1(it) + REAL(tab) * dist
             ENDIF
          ENDDO
          IF (ic1(it).EQ.0) dist1(it) = dist1(it) + dist
        ENDDO
      ENDDO
c...  Sort intesections:
      iccell = 0
      itcell = 0
      clcell = 0.0
      DO it =     osmnode(ind1)%tube_range(1), 
     .        MIN(osmnode(ind1)%tube_range(2),ntube)
        IF (tube(it)%type.EQ.GRD_BOUNDARY.OR.ic1(it).EQ.0) CYCLE
c        WRITE(0,*) 'PICKENS:',ir,ik1(ir),dist1(ir)
        IF     (it.EQ.itgive) THEN
          iccell(1) = ic1(it)   ! Current cell
          itcell(1) = it
        ENDIF
        IF (clcell(2).EQ.0.0.OR.dist1(it).LT.clcell(2)) THEN
          iccell(2) = ic1(it)
          itcell(2) = it
          clcell(2) = dist1(it) ! Cell farthest 'in'
        ENDIF
        IF (clcell(3).EQ.0.0.OR.dist1(it).GT.clcell(3)) THEN
          iccell(3) = ic1(it)
          itcell(3) = it
          clcell(3) = dist1(it) ! Cell farthest 'out'
        ENDIF
      ENDDO

c      WRITE(0,*) 'ICDATA:',iccell(1),iccell(2),iccell(3)
c      WRITE(0,*) '      :',itcell(1),itcell(2),itcell(3)
c      WRITE(0,*) '      :',cell(iccell(1))%ik,itcell(1)
c      WRITE(0,*) '      :',ind0,ind1
c      STOP 'sdfsd'
         
      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadUpstreamData(fname,fformat,mode,itube,coord,shift,
     .                            xcolumn,ycolumn,yval,s)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      CHARACTER, INTENT(IN)  :: fname*(*)
      INTEGER  , INTENT(IN)  :: fformat,itube,coord,mode
      REAL     , INTENT(IN)  :: shift,xcolumn,ycolumn,s
      REAL     , INTENT(OUT) :: yval

      INTEGER cind1,cind2
      LOGICAL osmGetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

      INTEGER    MAXTDAT     ,MAXCOLS   ,MAXNBUFFER_ARRAY    
      PARAMETER (MAXTDAT=1000,MAXCOLS=20,MAXNBUFFER_ARRAY=100)

      INTEGER   ndata,fp,i,j,n,ncolumns,xcol,ycol,ion,idum
      REAL      vdata(MAXTDAT,MAXCOLS),xval,frac
      CHARACTER     buffer*1024
      CHARACTER*256 buffer_array(MAXNBUFFER_ARRAY),buffer_temp

      TYPE(type_tube )              :: ass_tube
      TYPE(type_fluid), ALLOCATABLE :: ass_fluid(:,:)



c      WRITE(0,*) 'UPSTREAM: COORD2=',coord

c      WRITE(0,*) 'UPSTREAM: FILE='//TRIM(fname)//'<',mode

 
      IF (TRIM(fname).EQ.'reference_solution') THEN 
        ass_tube = tube(itube)

        cind1 = tube(itube)%cell_index(LO)
        cind2 = tube(itube)%cell_index(HI)          

        ALLOCATE(ass_fluid(cind2-cind1+1,nion))

        CALL InterpolateReferencePlasma
     .         (ass_tube,nion,ass_fluid,cell(cind1:cind2))

        ycol = NINT(ycolumn)

        ion = 1

        IF (mode.EQ.-1.OR.mode.EQ.-2) THEN 

          SELECTCASE (ycol)
            CASE (1)
              yval = ass_tube%jsat(-mode,ion) 
            CASE (4)
              yval = ass_tube%te  (-mode) 
            CASE (5)
              yval = ass_tube%ti  (-mode,ion) 
            CASE DEFAULT
              CALL ER('LoadUpstreamData','Unknown YCOLUMN',*99)
          ENDSELECT

        ELSE
c         Interpolate from the reference solution based on the distance
c         along the flux-tube, s:

          DO i = cind1, cind2-1
            IF (cell(i)%s.LT.s.AND.cell(i+1)%s.GT.s) EXIT
          ENDDO
          IF (i.EQ.cind2) 
     .      CALL ER('LoadUpstreamData','S not found',*99)

          j = i - cind1 + 1

          frac = (s - cell(i)%s) / (cell(i+1)%s - cell(i)%s)

          write(88,*) 'data dumper:',i,j,frac
          write(88,*) 'data dumper:',cind1,cind2
          write(88,*) '           :',s


          SELECTCASE (ycol)
            CASE (1)
              yval = (1.0 - frac) * ass_fluid(j  ,ion)%ne + 
     .                      frac  * ass_fluid(j+1,ion)%ne
            CASE (2)
              yval = (1.0 - frac) * ass_fluid(j  ,ion)%vi + 
     .                      frac  * ass_fluid(j+1,ion)%vi
            CASE (4)
              yval = (1.0 - frac) * ass_fluid(j  ,ion)%te + 
     .                      frac  * ass_fluid(j+1,ion)%te
            CASE (5)
              yval = (1.0 - frac) * ass_fluid(j  ,ion)%ti + 
     .                      frac  * ass_fluid(j+1,ion)%ti
            CASE DEFAULT
              CALL ER('LoadUpstreamData','Unknown YCOLUMN',*99)
          ENDSELECT

        ENDIF

        DEALLOCATE(ass_fluid)
        RETURN
      ENDIF

c...  Access data file:
      fp = 99
      OPEN(UNIT=fp,FILE=TRIM(fname),ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)

      DO WHILE (.TRUE.)
        READ(fp,*) buffer
        IF (buffer(1:1).NE.'*') EXIT
      ENDDO

      REWIND (fp)

      IF (buffer(1:1).EQ.'{') THEN 
        DO WHILE (osmGetLine(fp,buffer,WITH_TAG))
c          WRITE(0,*) 'BUFFER:',TRIM(buffer)
c...      Isolate tag string:
          DO i = 2, LEN_TRIM(buffer)
            IF (buffer(i:i).EQ.'}') EXIT
          ENDDO
          n = LEN_TRIM(buffer)
          SELECTCASE (buffer(2:i-1))
            CASE ('NUMBER OF COLUMNS')
              READ(buffer(i+1:n),*) ncolumns
            CASE ('DATA','DATA LIST')
              j = 0
              DO WHILE(osmGetLine(fp,buffer,NO_TAG))
                IF (buffer(1:1).EQ.'{') EXIT
                j = j + 1
                READ(buffer,*) vdata(j,1:ncolumns)
              ENDDO
              ndata = j
            CASE ('END')
              EXIT
            CASE DEFAULT
              CALL ER('LoadUpstreamData','Unknown tag',*99)
          ENDSELECT

        ENDDO
      ELSE

        ncolumns = -1
        ndata    =  0
        DO WHILE (.TRUE.)
          READ(fp,'(A1024)',END=10) buffer
c          write(0,*) 'buffer >'//buffer(1:20)//' <'
          IF (buffer(1:1).EQ.'*') CYCLE
          DO i = 1, MAXNBUFFER_ARRAY
            WRITE(buffer_array(i),'(256X)')
          ENDDO
          CALL SplitBuffer(buffer,buffer_array)
          idum = 0
          DO i = 1, MAXNBUFFER_ARRAY
            buffer_temp = buffer_array(i)
            IF (buffer_temp(1:1).NE.' ') idum = idum + 1
          ENDDO
c          WRITE(0,*) 'idum',idum
          IF (idum.EQ.0) CYCLE
          IF (ncolumns.EQ.-1) ncolumns = idum
          IF (idum.NE.ncolumns) THEN
            CALL ER('LoadUpstreamData','Number of columns in data '//
     .              'file changed',*99)
          ELSE
            ndata = ndata + 1
            DO i = 1, ncolumns
              READ(buffer_array(i),*) vdata(ndata,i)
            ENDDO
          ENDIF
        ENDDO
 
      ENDIF
10    CONTINUE
      CLOSE(fp)

c... 
      SELECTCASE (coord)
        CASE (3,5)
          xval = tube(itube)%psin - shift 
        CASE (4)
          xval = tube(itube)%rho - shift
        CASE DEFAULT
          CALL ER('LoadUpstreamData','Unknown COORD value',*99)
      ENDSELECT

      xcol = NINT(xcolumn)
      ycol = NINT(ycolumn)

      IF (xcol.LT.1.OR.xcol.GT.ncolumns) 
     .  CALL ER('LoadUpstreamData','x-data column ID invalid',*99)
      IF (ycol.LT.1.OR.ycol.GT.ncolumns.OR.xcol.EQ.ycol) 
     .  CALL ER('LoadUpstreamData','y-data column ID invalid',*99)

      IF     (xval.LT.vdata(1    ,xcol)) THEN
c        WRITE(0,*) 'WARNING:  INTERPOLATION FAILED, X-DATA BEYOND RANGE'
        yval = vdata(1,ycol)
      ELSEIF (xval.GT.vdata(ndata,xcol)) THEN
c        WRITE(0,*) 'WARNING:  INTERPOLATION FAILED, X-DATA BEYOND RANGE'
        yval = vdata(ndata,ycol)
      ELSE
        WRITE(88,*) 'FIT:',ndata
        WRITE(88,*) '   :',vdata(1:ndata,xcol)
        WRITE(88,*) '   :',vdata(1:ndata,ycol)

        CALL Fitter(ndata,vdata(1,xcol),vdata(1,ycol),
     .              1,xval,yval,'LINEAR')

        WRITE(88,*) '   :',xval,yval
        WRITE(88,*) '   : unshifted xval',xval + shift
      ENDIF

c      WRITE(88,*) 'XCOL,YCOL:',xcol,ycol
c      WRITE(88,*) 'SHIFT:',shift
c      WRITE(88,*) 'ITUBE,XVAL,YVAL:',itube,xval,yval

      RETURN
 98   WRITE(0,*) 'ERROR LoadUpstreamData: Data file not found'
 99   WRITE(0,*) '  FILE = ',TRIM(fname)
      WRITE(0,*) '  XCOL = ',xcol
      WRITE(0,*) '  YCOL = ',ycol
      WRITE(0,*) '  NCOL = ',ncolumns
      STOP
      END
c     
c ======================================================================
c
c Code from http://coding.derkeiler.com/Archive/Fortran/comp.lang.fortran/2005-01/0953.html
c
c
c      SUBROUTINE SetBit(word,position,value)
c      IMPLICIT none
c      INTEGER, INTENT(OUT) :: word
c      INTEGER, INTENT(IN)  :: position,value
c     
c      INTEGER :: bits(0:31)
c      
c      IF (value.NE.0.AND.value.NE.1) 
c     .  CALL ER('SetBit','Invalid VALUE for bit',*99)
c      IF (position.LT.0.OR.position.GT.31) 
c     .  CALL ER('SetBit','Invalid POSITION for bit',*99)
c
c      bits = IAND(ISHFT(word,[-31:0]),1)  ! unpack WORD
c
c      bits(position) = value
c
c      word = SUM(ISHFT(bits,[31:0],-1]))   ! pack bits into WORD
c
c      RETURN
c 99   WRITE(0,*) ' VALUE   =',value
c      WRITE(0,*) ' POSITION=',position
c      STOP
c      END
c     
c ======================================================================
c
c      LOGICAL FUNCTION GetBit(word,position)
c      IMPLICIT none
c      INTEGER, INTENT(IN) :: word,position
c     
c      INTEGER :: bits(0:31)
c      
c      IF (position.LT.0.OR.position.GT.31) 
c     .  CALL ER('GetBit','Invalid POSITION for bit',*99)
c
c      bits = IAND(ISHFT(word,[-31:0]),1)  ! unpack WORD
c
c      IF (bits(position).EQ.1) THEN
c        GetBit = .TRUE.
c      ELSE
c        GetBit = .FALSE.
c      ENDIF
c
c      RETURN
c 99   WRITE(0,*) ' POSITION=',position
c      STOP
c      END
c
c ======================================================================
c
      SUBROUTINE AssignNodeValues_New(itube,nnode,mnode,node,opt_tube)
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_sol28_targets
      USE mod_geometry
      USE mod_legacy
      IMPLICIT none

      INTEGER, INTENT(IN)  :: itube
      INTEGER, INTENT(OUT) :: nnode,mnode       
      TYPE(type_node), INTENT(OUT) :: node(*)
      TYPE(type_options_osm) :: opt_tube

      INTEGER CalcPoint
      REAL    GetJsat2,GetRelaxationFraction,GetCs2,CalcPressure

      INTEGER, PARAMETER :: MAXNNODES = 50

      INTEGER ic,ic1,ic2,it,it1,i0,i1,i2,i3,i4,ifit,index,mode,type,
     .        iobj,isrf,ivtx(2),nfit,icell(3),ion,coord,inode,iside,
     .        itarget,hold_ic,i5,check
      CHARACTER dummy*1024
      LOGICAL nc,vc,pc,tec,tic,density,tetarget,debug,link,intersection,
     .        first_pass,two_timer,default_message,node_valid,
     .        suppress_screen,rho_warning,
     .        ne_warning,te_warning,ti_warning
      REAL    te(0:6),ne(0:6),s(0:6),pe(0:6),ti(0:6),vb(0:6),
     .        frac,te0,te1,ti0,ti1,n0,n1,A,B,C,expon,te_cs,ti_cs,
     .        psin0,psin1,psin2,p0,p1,result(3),dist,radvel,rho0,
     .        prb1,tmp1,val,val0,val1,val2,p(0:5),v0,v1,v2,
     .        hold_tab,hold_tcd,ne_LO,ne_HI,cs,
     .        node_pe,node_v,node_ne,node_te,node_ti,smax,L,
     .        nustar 
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd,e1,e2,f1,f2,
     .        hold_c1,hold_c2,hold_d1,hold_d2,pts(3,4)


      INTEGER node_n,node_i(0:MAXNNODES)
      TYPE(type_node) :: node_s(0:MAXNNODES)

      TYPE(type_tube )              :: tube_tmp
      TYPE(type_fluid), ALLOCATABLE :: fluid_tmp(:,:)

      DATA default_message, rho_warning, ne_warning, te_warning, 
     .                                   ti_warning
     .     / .TRUE.       , .TRUE.     , .TRUE.    , .TRUE.    , 
     .                                   .TRUE.   /
      SAVE     

      suppress_screen = .FALSE.


      debug = .TRUE. 

      node_n = 0

      inode = 0

      ion = 1

      node_s%s  = 0.0
      node_s%ne = 0.0
      node_s%v  = 0.0
      node_s%pe = 0.0
      node_s%te = 0.0
      DO i1 = 1, MAXNNODES
        node_s(i1)%jsat(1:S28_MAXNION) = 0.0
        node_s(i1)%ti  (0:S28_MAXNION) = 0.0
      ENDDO

      frac = GetRelaxationFraction()

c...  Flag if ITUBE has been processed once already and was assigned a 
c     default symmetry point, since no proper point was found:
      two_timer = ibits(tube2(itube)%state,0,1).EQ.1  

c...  Better/cleaner to pass the tube to this routine, and not need to
c     use mod_sol28_locals...?
      it = itube

      write(logfp,*) '---> here in assignnodevalues',tube(it)%n


      DO i1 = 2, osmnnode
        i0 = i1 - 1

        node_valid = .TRUE.

        IF (debug) THEN
          WRITE(logfp,*) 'INTER:',i0,i1,osmnode(i0)%type,
     .                            osmnode(i1)%type
          WRITE(logfp,*) 'PAR_MODE:',osmnode(i0:i1)%par_mode
        ENDIF

        IF (osmnode(i0)%type.NE.osmnode(i1)%type.OR.
     .      osmnode(i1)%type.EQ.0.0.OR.
     .      osmnode(i1)%type.EQ.0.0) THEN
          first_pass = .TRUE.
          CYCLE
        ENDIF

        index = 0

        ne = 0.0  ! Necessary, also appears below...
        vb = 0.0
        pe = 0.0
        te = 0.0
        ti = 0.0
c...    

        IF (debug) THEN
          WRITE(logfp,*) 'INTER GO:',it,osmnode(i1)%tube_range(1:2)
        ENDIF

c...    Do not apply data if IT is outside specified range:
        IF ((it.LT.osmnode(i1)%tube_range(1).OR.
     .       it.GT.osmnode(i1)%tube_range(2)).AND.
     .      osmnode(i1)%type.NE.3.0) CYCLE

c...    Need to keep track of iteration currently being assigned, so
c       that scans in upstream with step can be done serially, rather than
c       having everything on the same line...

c...    Check that rings from different grid regions are not in the same group
c       of rings:
        DO i2 =     osmnode(i1)%tube_range(1), 
     .          MIN(osmnode(i1)%tube_range(2),ntube)-1
          IF (i2.GT.ntube-1) EXIT          
          IF (tube(i2)%type.NE.tube(i2+1)%type) THEN
            IF (logop.GT.0.AND.tube(i2)%type.NE.GRD_CORE) THEN
              WRITE(logfp,*)
              WRITE(logfp,*) '-------------------------------------'
              WRITE(logfp,*) ' THAT FUNNY THING ABOUT MIXED REG.!? '
              WRITE(logfp,*) '--------------------------- ---------'
              WRITE(logfp,*)
            ENDIF          
          ENDIF      
        ENDDO

c MODE                          P1           P2
c   1 - power law               coordinate   index
c   2 - exponential v0-v2       coordinate   index
c   3 - exponential to infinity
c   4 - from probe data         coordinate   probe number
c   5 - parameter fits
c   6 - core + pedestal + SOL automated fits
c   7 - exponential decay for velocity, temperature, decay based on v_perp and L for density
c

c   coord = 1 - linear on line segment
c         = 2 - RHO over the range of applicability
c         = 3 - PSIn over range of applicability (like coord=2) 
c         = 4 - RHO
c         = 5 - PSIn (raw)
c         = 6 - linear on line segment, but from tube link to infinity
c         = 7 - psin, from tube link to infinity

        index = NINT(osmnode(i1)%type)
        mode  = osmnode(i1)%rad_mode
        coord = osmnode(i1)%rad_coord
        expon = osmnode(i1)%rad_exp

        IF (debug) THEN
          WRITE(logfp,'(6A)') '  I1',' INDEX','  MODE',' COORD',' EXPON'
          WRITE(logfp,'(6A)') '  --',' -----','  ----',' -----',' -----'
          WRITE(logfp,'(I4,3I6,F6.1)') i1,index,mode,coord,expon
        ENDIF

c...    Decide if specified upstream data is density or pressure:
        density = .TRUE.   ! *** I DON'T LIKE THIS SOLUTION ***
c        IF (index.EQ.3.AND.
c     .      (tube(it)%type.EQ.GRD_SOL.AND..FALSE..OR.
c     .       tube(it)%type.EQ.GRD_PFZ.AND..FALSE.))
c     .    density = .FALSE.

        IF (mode.EQ.7.AND..NOT.two_timer) THEN
          suppress_screen = .TRUE.
          a1 = 0.0
          a2 = 0.0
          b1 = 0.0
          b2 = 0.0
        ELSE
          a1 = DBLE(osmnode(i1-1)%rad_x)
          a2 = DBLE(osmnode(i1-1)%rad_y)
          b1 = DBLE(osmnode(i1  )%rad_x)
          b2 = DBLE(osmnode(i1  )%rad_y)
        ENDIF

        intersection = .FALSE.

        IF (osmnode(i1)%type.EQ.3.0) THEN
c...      Target nodes:
          IF (osmnode(i1)%tube_range(1).NE.osmnode(i1)%tube_range(2))
     .      CALL ER('AssignNode_New','Single target position must '//
     .              'be specified',*99)
          DO i4 = 1, ntarget
            IF (target(i4)%location.EQ.ABS(osmnode(i1)%tube_range(1))) 
     .        EXIT
          ENDDO
          IF (i4.LT.ntarget+1) THEN         
            IF (osmnode(i1)%tube_range(1).LT.0) THEN
c             All tubes on grid assigned:
              i5 = 0  
            ELSE
c             Only select tubes identified as being part of this target data block:
              DO i5 = 1, target(i4)%nlist
                IF (target(i4)%ilist(i5).EQ.it) EXIT
              ENDDO
            ENDIF
            IF (i5.LT.target(i4)%nlist+1) THEN
              intersection = .TRUE.  
              hold_c1  = 0.0D0
              hold_c2  = 0.0D0
              hold_d1  = 0.0D0
              hold_d2  = 0.0D0
              hold_tab = 0.0D0
              IF (debug) THEN
                WRITE(logfp,*) 'Target assignment identified:',i4
                WRITE(logfp,*) ' SPEC L=',osmnode(i1)%tube_range(1:2),i1
                WRITE(logfp,*) ' LOC   =',target(i4)%location
                WRITE(logfp,*) ' TAG   =',TRIM(target(i4)%tag)
                WRITE(logfp,*) ' POS   =',target(i4)%position
              ENDIF
              IF (target(i4)%position.EQ.LO) THEN
                hold_ic  = tube(it)%cell_index(LO)
                hold_tcd = 0.0D0
                osmnode(i0:i1)%par_mode = -1
              ELSE
                hold_ic  = tube(it)%cell_index(HI)
                hold_tcd = 1.0D0
                osmnode(i0:i1)%par_mode = -2
              ENDIF
            ENDIF
          ENDIF
c       ----------------------------------------------------------------  
        ELSEIF (osmnode(i1)%type.EQ.2.1.AND.             ! Semi-automated symmetry point specification
     .          two_timer) THEN                          ! WHAT IS THIS!? -SL, 17/10/2013                    
          intersection = .FALSE.
          hold_c1  = 0.0D0
          hold_c2  = 0.0D0
          hold_d1  = 0.0D0
          hold_d2  = 0.0D0
          hold_tab = 0.0D0
c...      Check if there's an interesection (optional):
          hold_ic = -1
          DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)
            iobj = ic
            isrf = ABS(obj(iobj)%iside(1))          ! *** Use GetVertex *** 
            ivtx(1:2) = srf(isrf)%ivtx(1:2)
            c1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
            c2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2)))
            isrf = ABS(obj(iobj)%iside(3))
            ivtx(1:2) = srf(isrf)%ivtx(1:2)
            d1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
            d2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2)))
            CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
            IF (tab.GE.0.0D0.AND.tab.LT.1.0D0.AND.
     .          tcd.GE.0.0D0.AND.tcd.LT.1.0D0) THEN
              intersection = .TRUE.
              hold_ic  = ic 
              IF (debug) THEN
                WRITE(logfp,*) '  intersection'
                WRITE(logfp,*) '  tcd ',tcd
                WRITE(logfp,*) '  c1,2',c1,c2
                WRITE(logfp,*) '  d1,2',d1,d2
              ENDIF
              hold_c1  = c1 
              hold_c2  = c2
              hold_d1  = d1
              hold_d2  = d2
              hold_tab = tab
              hold_tcd = tcd
              EXIT
            ENDIF
          ENDDO
c...      No interesection found, so take the middle of the tube (in s):
          IF (hold_ic.EQ.-1) THEN           
            IF (debug)
     .        write(logfp,*) 'debug: trying to find...'
            DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)            
              IF (cell(ic)%sbnd(1).LE.0.5*tube(it)%smax.AND.
     .            cell(ic)%sbnd(2).GE.0.5*tube(it)%smax) THEN
                IF (debug)
     .            write(logfp,*) 'debug: made up interesecton',ic,
     .                           tube(it)%cell_index(LO:HI)
                hold_ic = ic
                EXIT
              ENDIF
            ENDDO
            IF (ic.EQ.tube(it)%cell_index(HI)+1)
     .        CALL ER('AssignNodeValues_New','Two-timing symmetry '//
     .                'node cell not identified',*99)
          ENDIF
c       ----------------------------------------------------------------  
        ELSE
          DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)
c...        Assumed 1:1 mapping between grid and data:
            iobj = ic
            isrf = ABS(obj(iobj)%iside(1))          ! *** Use GetVertex *** 
            ivtx(1:2) = srf(isrf)%ivtx(1:2)
            c1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
            c2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2)))
            isrf = ABS(obj(iobj)%iside(3))
            ivtx(1:2) = srf(isrf)%ivtx(1:2)
            d1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
            d2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2)))
            CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
            IF ((tab.GE.0.0D0.AND.tab.LT.1.0D0.AND.
     .           tcd.GE.0.0D0.AND.tcd.LT.1.0D0).OR.
     .          (ic.EQ.tube(it)%cell_index(LO).AND.
     .           (osmnode(i0)%par_mode.EQ.-1.OR.
     .            osmnode(i0)%par_mode.EQ.-2))) THEN
              intersection = .TRUE.
              hold_ic  = ic
              hold_c1  = c1  ! why not just trigger an exit here and avoid the whole "hold_" business?
              hold_c2  = c2  ! (also find this above for compatibility, but with an EXIT)
              hold_d1  = d1
              hold_d2  = d2
              hold_tab = tab
              hold_tcd = tcd
            ENDIF
          ENDDO
        ENDIF

        IF (.NOT.intersection.AND.
     .      .NOT.(osmnode(i1)%type.EQ.2.1.AND.two_timer)) CYCLE

c...    An intersection between the line segment and the ring has been found:

        IF (.NOT.first_pass) THEN
          WRITE(0,*) 'ERROR: interpolation node block intersects '//
     .               'the flux-tube more than once, are you mad'
          WRITE(0,*) '  ITUBE = ',itube
          WRITE(0,*) '  I0,1  = ',i0,i1
          STOP
        ENDIF

        s  = 0.0
        ne = 0.0
        vb = 0.0
        pe = 0.0
        te = 0.0
        ti = 0.0
        link = .FALSE.

        ic  = hold_ic 
        c1  = hold_c1
        c2  = hold_c2
        d1  = hold_d1
        d2  = hold_d2
        tab = hold_tab
        tcd = hold_tcd

        IF (debug) WRITE(logfp,*) 
     .    'INTERSECTION:',intersection,i0,i1,itube,ic

c...    These are here in case psin0 gets picked up when linking exponential
c       decay data to neighbouring ring:
        IF (tube(it)%type.EQ.GRD_PFZ) THEN
          psin0 = -RHI
          psin1 =  RHI
          rho0  =  RHI
        ELSE
          psin0 =  RHI
          psin1 =  RHI
          rho0  =  RHI
        ENDIF

        IF (index.LT.1.OR.index.GT.3) 
     .    CALL ER('AssignNodeValues_2','Invalid parameter '//
     .            'index',*99)

        inode = 1 

        IF (intersection) THEN
          s(inode) = cell(ic)%sbnd(1) + SNGL(tcd) * cell(ic)%ds
        ELSE
          s(inode) = cell(ic)%s
        ENDIF

c...    Find data boundary values -- NEEDS WORK!:
        i2 = i0
        i3 = i1

        IF     (mode.EQ.1.OR.mode.EQ.2.OR.mode.EQ.3.OR.
     .          mode.EQ.5.OR.mode.EQ.6.OR.mode.EQ.7) THEN
c...      Interpolation boundary values provided in the input file:

c         *CRAP!*
          n0  = osmnode(i2)%ne  ! Needs work to allow relaxation 
          n1  = osmnode(i3)%ne
          v0  = osmnode(i2)%v
          v1  = osmnode(i3)%v
          p0  = osmnode(i2)%pe
          p1  = osmnode(i3)%pe
          te0 = osmnode(i2)%te
          te1 = osmnode(i3)%te
          ti0 = osmnode(i2)%ti(1)
          ti1 = osmnode(i3)%ti(1)

          IF (debug) THEN
            WRITE(logfp,*) 's,smax:',s(inode),tube(it)%smax
            WRITE(logfp,*) 'N0,1  :',n0,n1
            WRITE(logfp,*) 'V0,1  :',v0,v1
            WRITE(logfp,*) 'P0,1  :',p0,p1
            WRITE(logfp,*) 'Te0,1 :',te0,te1
            WRITE(logfp,*) 'Ti0,1 :',ti0,ti1
          ENDIF

          IF (osmnode(i2)%ne.EQ.-99.0.OR.
     .        osmnode(i2)%v .EQ.-99.0.OR.
     .        osmnode(i2)%pe.EQ.-99.0.OR.
     .        osmnode(i2)%te.EQ.-99.0) THEN
c...        Linking to another plasma region where the solution has
c           already (!) been calculated: 
            link = .TRUE.

            IF (osmnode(i1)%type.EQ.2.1.AND.two_timer) THEN
              icell = hold_ic
            ELSE
              CALL FindCell_New(i2,i3,it,icell) ! ,e1,e2) -bug, SL, 08/10/2010
            ENDIF

c...        Assumptions: 1:1 mapping between cells and objects, the
c           objects are 4 sided and have the standard DIVIMP indexing.  
c           Also, cells in the core/SOL will always reference tubes
c           that have a lower tube index via side 1-4 and PFZ cells 
c           reference tubes with a higher index through side 2-3.
            iobj = icell(2)
            IF (tube(it)%type.EQ.GRD_PFZ) THEN
              iside = 2
            ELSE
              iside = 4
            ENDIF                 
            ic1 = obj(iobj)%omap(iside)
            IF (debug) THEN
              WRITE(logfp,*) ' MAP   ',iobj,iside,ic1,nobj,ncell
              WRITE(logfp,*) ' MAP   ',tube(itube)%cell_index(LO:HI)
              WRITE(logfp,*) ' MAP   ',it
              WRITE(logfp,*) ' MAP   ',icell
              WRITE(logfp,*) ' MAP O ',
     .           obj(tube(it)%cell_index(LO):
     .               tube(it)%cell_index(HI))%omap(4)
            ENDIF
c...        Now have to search and see which tube this cell is in:
            IF (ic1.LE.0) THEN
c            IF (ic1.EQ.0) THEN
c...          The connection map has failed, perhaps due to generalised
c             DIVIMP grids, so map the cell on the fly:
              CALL DynamicMap(iobj,iside,ic1,it1)
            ELSE
              DO it1 = 1, ntube   ! *** REPLACE WITH GETTUBE CALL ***
c                IF (debug) WRITE(logfp,*) 
c     .           'TUBES:',it1,tube(it1)%cell_index(LO:HI)
                IF (tube(it1)%cell_index(LO).LE.ic1.AND.
     .              tube(it1)%cell_index(HI).GE.ic1) EXIT
              ENDDO
            ENDIF
            IF (it1.EQ.-1.OR.it1.EQ.ntube+1) 
     .        CALL ER('AssignNodeValues_New','Tube not identified',*99)

            IF (two_timer.AND.ibits(tube2(it1)%state,0,1).EQ.1) THEN
              STOP 'TRYING TO REFERENCE A DEFAULT TUBE'
            ENDIF            

            IF (osmnode(i2)%ne.EQ.-99.0) THEN
              IF (density) THEN
                n0 = fluid(ic1,1)%ne
              ELSE
                n0 = 2.0 * fluid(ic1,1)%te * fluid(ic1,1)%ne  ! Add Ti and M? 
              ENDIF
            ENDIF
            IF (osmnode(i2)%v    .EQ.-99.0) v0 = fluid(ic1,1)%vi + 0.003
            IF (osmnode(i2)%pe.EQ.-99.0)   ! No v|| contribution!
     .        p0 = fluid(ic1,1)%ne * 
     .             (fluid(ic1,1)%te + fluid(ic1,1)%ti)
            IF (osmnode(i2)%te   .EQ.-99.0) te0 = fluid(ic1,1)%te
            IF (osmnode(i2)%ti(1).EQ.-99.0) ti0 = fluid(ic1,1)%ti

c...        Base second radial interpolation value on the first value:
            IF (osmnode(i3)%ne   .LT.  0.0.AND.
     .          osmnode(i3)%ne   .NE.-77.0) n1 = -osmnode(i3)%ne   * n0
            IF (osmnode(i3)%pe   .LT.  0.0.AND.
     .          osmnode(i3)%pe   .NE.-77.0.AND.
     .          osmnode(i3)%pe   .NE.-88.0) p1 = -osmnode(i3)%pe   * p0
            IF (osmnode(i3)%te   .LT.  0.0.AND.
     .          osmnode(i3)%te   .NE.-77.0) te1= -osmnode(i3)%te   * te0
            IF (osmnode(i3)%ti(1).LT.  0.0.AND.
     .          osmnode(i3)%ti(1).NE.-77.0) ti1= -osmnode(i3)%ti(1)* ti0

            IF (coord.EQ.3.OR.coord.EQ.7) psin0 = tube(it1)%psin
            IF (coord.EQ.8)               rho0  = tube(it1)%rho  

            IF (debug) THEN
              WRITE(logfp,*) 'NODE LINK:',ic1,it1
              WRITE(logfp,*) '    IK,IR:',cell(ic1)%ik,cell(ic1)%ir
              WRITE(logfp,*) '    psin0:',psin0
              WRITE(logfp,*) '     rho0:',rho0
              WRITE(logfp,*) '       ne:',n0  
              WRITE(logfp,*) '       v :',v0  
              WRITE(logfp,*) '       pe:',p0  
              WRITE(logfp,*) '      Te0:',te0  
              WRITE(logfp,*) '      Ti0:',ti0  
            ENDIF

c...        Check if the tube being linked to has already been processed  -- MOD moved this here from just below the IF block ---
c           this iteration:                                                  SL, 26/01/2011
            IF (ibits(tube2(it1)%state,1,1).EQ.0) THEN  
              write(0    ,*) 'linking a tube that is not defined'//
     .                       ', bad',itube
              write(logfp,*) 'linking a tube that is not defined'//
     .                       ', bad',itube
              tube2(itube)%state = ibset(tube2(itube)%state,2)
              node_valid = .FALSE.
            ENDIF

          ENDIF


c...      Make sure that te0,1 are positive:
          IF (tetarget) THEN
            te0 = ABS(te0)
            te1 = ABS(te1)
          ENDIF

          IF     (coord.EQ.1) THEN
c...        Linear along the line segment, nice and simple:
            val0 = 0.0
            val1 = 1.0
            val = SNGL(tab)
          ELSEIF (coord.EQ.2) THEN
            val0 = tube(osmnode(i1)%tube_range(1))%rho
            val1 = tube(osmnode(i1)%tube_range(2))%rho
            val  = tube(itube)%rho
          ELSEIF (coord.EQ.4) THEN
c...        Well, some funny business here since I can't imagine why I took the absolute value 
c           of RHO, which is causing a discrepancy between i-ref-0001d and i-ref-0007d, where the
c           latter is a re-run of the former with iter_steven_3 (and the former was iter_steven_2).
c           This change is also in iter_steven_2 however, in contrast to when i-ref-0001d was run,
c           which had VAL<0 in the core, as did i-ref-0009d, so this must be a relatively recent
c           change.  Very odd.
            val = tube(it)%rho
c            val = ABS(tube(it)%rho)
            IF (val.EQ.0.0) 
     .        WRITE(0,*) 'WARNING AssignNodeValues: RHO=0.0 for '//
     .                   'tube',it
            IF (rho_warning) THEN 
              rho_warning = .FALSE.
              WRITE(logfp,*) 
              WRITE(logfp,*) '-----------------------------------------'
              WRITE(logfp,*) '   ABS() REMOVED FROM RHO ASSIGNMENT!    '
              WRITE(logfp,*) '-----------------------------------------'
              WRITE(logfp,*) 
            ENDIF
          ELSEIF (coord.EQ.5) THEN
c...        Just PSIn:
            psin0 = tube(it)%psin
            val   = psin0
            WRITE(logfp,*) ' psin 5:',psin0
          ELSEIF (coord.EQ.6) THEN
c...        Linear along the line segment, starting at the link:
            IF (.NOT.link) 
     .        CALL ER('AssignNodeValues_New','COORD=6 but no '//
     .                'link to calculated plasma',*99)
            val0 = 0.0
            val1 = 1.0
            val  = 1.0E+6
            IF (debug) WRITE(logfp,*) '  intersection',intersection
            IF (intersection) THEN
              f1 = c1 + tcd * (d1 - c1)  ! Point where the current focus tube intersects
              f2 = c2 + tcd * (d2 - c2)  ! the interpolation line segment
              IF (debug) THEN
                WRITE(logfp,*) '  tcd ',tcd
                WRITE(logfp,*) '  c1,2',c1,c2
                WRITE(logfp,*) '  d1,2',d1,d2
                WRITE(logfp,*) '  f1,2',f1,f2
              ENDIF
              DO ic2 = ic1-1, ic1+1
c...            Assuming a 1:1 mapping between grid and data:
                iobj = ic2                                        ! *** Replace with GETVERTEX! ***
                isrf = ABS(obj(iobj)%iside(1))
                ivtx(1:2) = srf(isrf)%ivtx(1:2)
                c1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
                c2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2)))
                isrf = ABS(obj(iobj)%iside(3))
                ivtx(1:2) = srf(isrf)%ivtx(1:2)
                d1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
                d2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2)))
                CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
                IF (tab.GE.0.0D0.AND.tab.LT.1.0D0.AND.
     .              tcd.GE.0.0D0.AND.tcd.LT.1.0D0) THEN
                  e1 = a1 + tab * (b1 - a1)  ! Point where the linked tube intesects 
                  e2 = a2 + tab * (b2 - a2)  ! the interpoloation line segment
                  val = SNGL(DSQRT((e1-f1)**2+(e2-f2)**2))
                  IF (debug) THEN
                    WRITE(logfp,*) '-> it1,ic2,range_it1',it1,ic2,
     .                             tube(it1)%cell_index(LO:HI)
                    WRITE(logfp,*) '  e1,2',e1,e2
                    WRITE(logfp,*) '  f1,2',f1,f2
                    WRITE(logfp,*) '  val ',val
                  ENDIF
                  EXIT
                ENDIF
              ENDDO
            ELSE
              f1 = cell(ic)%cencar(1)  ! Centre of the focus cell on the current tube
              f2 = cell(ic)%cencar(2)  
              IF (.TRUE.) THEN 
                DO ic1 = tube(it1)%cell_index(LO),
     .                   tube(it1)%cell_index(HI)
                  iobj = ic1 
c                 Get the centreline of the cell:
                  CALL GetVertex(iobj,1,pts(1,1),pts(2,1))
                  CALL GetVertex(iobj,2,pts(1,2),pts(2,2))
                  c1 = 0.5D0 * (pts(1,1) + pts(1,2))
                  c2 = 0.5D0 * (pts(2,1) + pts(2,2))
                  CALL GetVertex(iobj,3,pts(1,1),pts(2,1))
                  CALL GetVertex(iobj,4,pts(1,2),pts(2,2))
                  d1 = 0.5D0 * (pts(1,1) + pts(1,2))
                  d2 = 0.5D0 * (pts(2,1) + pts(2,2))
                  c1 = 0.5D0 * (c1 + d1)
                  c2 = 0.5D0 * (c2 + d2)
                  dist = SNGL(DSQRT( (f1-c1)**2 + (f2-c2)**2 ))
                  IF (debug) THEN
                    WRITE(logfp,*) '  check ',check,tcd
                    WRITE(logfp,*) '        ',c1,c2
                    WRITE(logfp,*) '        ',f1,f2
                    WRITE(logfp,*) '        ',dist
                  ENDIF          
                  IF (dist.LT.val) THEN
                    IF (debug) THEN
                      WRITE(logfp,*) '  val update',dist
                    ENDIF
                    val = dist
                  ENDIF
                ENDDO               
              ELSE
c               ORIGINAL METHOD, BUT NOT FOOLPROOF
c...            Scan along the link tube and find the shortest perpendicular 
c               distance to the focus cell:
                DO ic1 = tube(it1)%cell_index(LO),
     .                   tube(it1)%cell_index(HI)
                  iobj = ic1 
c                 Get the centreline of the cell:
                  CALL GetVertex(iobj,1,pts(1,1),pts(2,1))
                  CALL GetVertex(iobj,2,pts(1,2),pts(2,2))
                  c1 = 0.5D0 * (pts(1,1) + pts(1,2))
                  c2 = 0.5D0 * (pts(2,1) + pts(2,2))
                  CALL GetVertex(iobj,3,pts(1,1),pts(2,1))
                  CALL GetVertex(iobj,4,pts(1,2),pts(2,2))
                  d1 = 0.5D0 * (pts(1,1) + pts(1,2))
                  d2 = 0.5D0 * (pts(2,1) + pts(2,2))
                  check = CalcPoint(c1,c2,d1,d2,f1,f2,tcd)
                  IF (debug) THEN
                    WRITE(logfp,*) '  check ',check,tcd
                    WRITE(logfp,*) '        ',c1,c2
                    WRITE(logfp,*) '        ',d1,d2
                    WRITE(logfp,*) '        ',f1,f2
                  ENDIF          
                  IF (check.EQ.2) THEN
                    e1 = c1 + tcd * (d1 - c1)  
                    e2 = c2 + tcd * (d2 - c2)  
                    dist = SNGL(DSQRT( (e1-f1)**2 + (e2-f2)**2 ))
                    WRITE(logfp,*) '  val check ',val,dist
                    WRITE(logfp,*) '            ',e1,e2
                    IF (dist.LT.val) THEN
                      IF (debug) THEN
                        WRITE(logfp,*) '  val update',dist
                      ENDIF
                      val = dist
                    ENDIF
                  ENDIF
                ENDDO
c                WRITE(0,*) 'itube,it1=',itube,it1
c                STOP 'sdfdsf'
              ENDIF
            ENDIF
            IF (val.EQ.1.0E+6) 
     .        CALL ER('AssignNodeValues_New','Linear link '//
     .                'reference not found',*99)

            IF (val.EQ.1.0E+6) 
     .        CALL ER('AssignNodeValues_New','Linear link '//
     .                'reference not found',*99)
            IF (debug) WRITE(logfp,*) '6: VAL=',val
          ELSEIF (coord.EQ.3.OR.coord.EQ.7) THEN
            val0 = 0.0  ! Necessary?
            val1 = 0.0  
            val  = ABS(tube(it)%psin - psin0)
            WRITE(logfp,*) ' psin 7:',psin0,tube(it)%psin,val
          ELSEIF (coord.EQ.8) THEN
            val0 = 0.0  ! Necessary?
            val1 = 0.0  
            val  = ABS(tube(it)%rho - rho0)
            WRITE(logfp,*) ' rho 8:',rho0,tube(it)%rho,val
          ELSE
c...         *** THIS IS REALLY OLD CODE I THINK -- EFFECTIVELY REPLACED ABOUVE BY FINDCELL_NEW... CHECK...
            STOP 'ASSIGNNODEVALUES: CODE IS OBSOLETE'
c           * big snip * -SL, 17/10/2013
          ENDIF
        ELSEIF (mode.EQ.4) THEN
c...      Load probe data, dummy values here:
          n0  = osmnode(i2)%ne
          p0  = osmnode(i2)%pe
          v0  = osmnode(i2)%v 
          te0 = osmnode(i2)%te
          ti0 = osmnode(i2)%ti(1)
          n1  = n0
          p1  = p0
          v1  = v0
          te1 = te0
          ti1 = ti0
        ELSE
          CALL ER('AssignNodeValues_2','Invalid MODE A',*99)   
        ENDIF

        IF (debug) WRITE(logfp,*) 'VAL:',val0,val1,val

c...    Check if quantities should be assigned:
        nc  = .TRUE.
        vc  = .TRUE.
        pc  = .TRUE.
        tec = .TRUE.
        tic = .TRUE.
        IF (n0 .EQ.0.0.OR.n1 .EQ.0.0) nc  = .FALSE.
        IF (v0 .EQ.0.0.OR.v1 .EQ.0.0) vc  = .FALSE.
        IF (p0 .EQ.0.0.OR.p1 .EQ.0.0) pc  = .FALSE.
        IF (te0.EQ.0.0.OR.te1.EQ.0.0) tec = .FALSE.
        IF (ti0.EQ.0.0.OR.ti1.EQ.0.0) tic = .FALSE.

        IF (debug) THEN
          WRITE(logfp,*) 'again'
          WRITE(logfp,*) 'N0,1 :',n0,n1
          WRITE(logfp,*) 'V0,1 :',v0,v1
          WRITE(logfp,*) 'P0,1 :',p0,p1
          WRITE(logfp,*) 'Te0,1:',te0,te1
          WRITE(logfp,*) 'Ti0,1:',ti0,ti1
        ENDIF

        SELECTCASE (mode)
c         ----------------------------------------------------------
          CASE (1)
c...        Power law between v1 and v2:
            val2 = (val - val0) / (val1 - val0)
            IF (nc ) ne(inode) = n0  + val2**expon * (n1  - n0 )
            IF (vc ) vb(inode) = v0  + val2**expon * (v1  - v0 )
            IF (pc ) pe(inode) = p0  + val2**expon * (p1  - p0 )
            IF (tec) te(inode) = te0 + val2**expon * (te1 - te0)
            IF (tic) ti(inode) = ti0 + val2**expon * (ti1 - ti0)
c         ----------------------------------------------------------
          CASE (2)
c...        Exponential decay between v1 and v2:
            C = expon  
            IF (nc) THEN
              A = (n1 - n0) / (EXP(-val1 / C) - 1.0)
              B = n0 - A
              ne(inode) = A * EXP(-val / C) + B
            ENDIF
            IF (vc) THEN
              A = (v1 - v0) / (EXP(-val1 / C) - 1.0)
              B = v0 - A
              vb(inode) = A * EXP(-val / C) + B
            ENDIF
            IF (pc) THEN
              A = (p1 - p0) / (EXP(-val1 / C) - 1.0)
              B = p0 - A
              pe(inode) = A * EXP(-val / C) + B
            ENDIF
            IF (tec) THEN
              A = (te1 - te0) / (EXP(-val1 / C) - 1.0)
              B = te0 - A
              te(inode) = A * EXP(-val / C) + B
            ENDIF
            IF (tic) THEN
              A = (ti1 - ti0) / (EXP(-val1 / C) - 1.0)
              B = ti0 - A
              ti(inode) = A * EXP(-val / C) + B
            ENDIF
c         ----------------------------------------------------------
          CASE (3)
c...        Exponential decay toward some value:
            IF (debug) WRITE(logfp,*) 'MODE=3: ',te0,te1,val
            C = expon
            A = n0 - n1
            B = n1
            IF (nc) ne(inode) = A * EXP(-val / C) + B
            C = expon
            A = v0 - v1
            B = v1
            IF (vc) vb(inode) = A * EXP(-val / C) + B
            A = p0 - p1
            B = p1
            IF (pc) pe(inode) = A * EXP(-val / C) + B
            WRITE(logfp,*) 'pc:',pc
            WRITE(logfp,*) '  :',p0,p1
            WRITE(logfp,*) '  :',A,B,val,C
            A = te0 - te1
            B = te1 
            IF (tec) te(inode) = A * EXP(-val / C) + B
            WRITE(logfp,*) 'tec:',tec
            WRITE(logfp,*) '   :',A,B,val,C
            A = ti0 - ti1
            B = ti1 
            IF (tic) ti(inode) = A * EXP(-val / C) + B
c         ----------------------------------------------------------
          CASE (4)
c...        Load probe data from ASCII file:
            WRITE(logfp,*) 'DATA FILE '//TRIM(osmnode(i1)%file_name)
            WRITE(logfp,*) 'COORD     ',coord

            IF (nc) THEN
              CALL LoadUpstreamData(osmnode(i1)%file_name,
     .               osmnode(i1)%file_format,
     .               osmnode(i1)%par_mode,
     .               itube,coord,osmnode(i1)%file_shift,
     .               expon,osmnode(i2)%ne,tmp1,s(inode))
              ne(inode) = tmp1 * osmnode(i1)%file_scale_ne
            ENDIF
            IF (vc) THEN
              CALL LoadUpstreamData(osmnode(i1)%file_name,
     .               osmnode(i1)%file_format,
     .               osmnode(i1)%par_mode,
     .               itube,coord,osmnode(i1)%file_shift,
     .               expon,osmnode(i2)%v ,tmp1,s(inode)) 
              vb(inode) = tmp1 * osmnode(i1)%file_scale_M
              WRITE(logfp,*) 'VB  B:',vb(inode)
            ENDIF
            IF (pc) THEN
              CALL LoadUpstreamData(osmnode(i1)%file_name,
     .               osmnode(i1)%file_format,
     .               osmnode(i1)%par_mode,
     .               itube,coord,osmnode(i1)%file_shift,
     .               expon,osmnode(i2)%pe,tmp1,s(inode)) 
              pe(inode) = tmp1 * osmnode(i1)%file_scale_pe
              WRITE(logfp,*) 'PE  B:',pe(inode)
            ENDIF
            IF (tec) THEN
              CALL LoadUpstreamData(osmnode(i1)%file_name,
     .               osmnode(i1)%file_format,
     .               osmnode(i1)%par_mode,
     .               itube,coord,osmnode(i1)%file_shift,
     .               expon,osmnode(i2)%te,tmp1,s(inode)) 
              te(inode) = tmp1 * osmnode(i1)%file_scale_te
            ENDIF
            IF (tic) THEN
              CALL LoadUpstreamData(osmnode(i1)%file_name,   ! Necessary to call LoadUpstreamData so many times?
     .               osmnode(i1)%file_format,
     .               osmnode(i1)%par_mode,
     .               itube,coord,osmnode(i1)%file_shift,
     .               expon,osmnode(i2)%ti(1),tmp1,s(inode)) 
              ti(inode) = tmp1 * osmnode(i1)%file_scale_ti
              WRITE(logfp,*) 'TI B:',ti(inode)
            ENDIF
c         ----------------------------------------------------------
          CASE (5)
c...        Interpolations from polynomial and exponential fitting parameters 
c           that are listed in the input file:
 
c...        Count how many fit data lines are to be processed in this group:
            nfit = 0
            DO i4 = i1+1, osmnnode
              IF (osmnode(i4)%type.NE.0.0) EXIT  
              nfit = nfit + 1
            ENDDO
c            WRITE(logfp,*) 'NFIT:',nfit

            DO i4 = 1, nfit
              ifit = i1 + i4

c              WRITE(logfp,*) 'FIT:',osmnode(ifit)%fit_type
c              WRITE(logfp,*) '   :',osmnode(ifit)%fit_psin(1)
c              WRITE(logfp,*) '   :',osmnode(ifit)%fit_psin(2)
c              WRITE(logfp,*) '   :',osmnode(ifit)%fit_shift
c              WRITE(logfp,*) '   :',osmnode(ifit)%fit_quantity
c              WRITE(logfp,*) '   :',osmnode(ifit)%fit_p(1:6)

              IF     (osmnode(ifit)%fit_type.EQ. 0.0) THEN 
              ELSEIF (osmnode(ifit)%fit_type.EQ.-1.0) THEN 
                psin1 = psin0 - osmnode(ifit)%fit_psin(1)
              ELSEIF (psin1.GE.osmnode(ifit)%fit_psin(1).AND.
     .                psin1.LE.osmnode(ifit)%fit_psin(2)) THEN
c...            Only need this for the exponential fit coefficients
c               returned in IDL (ts.pro at the moment):
                IF (osmnode(ifit)%fit_type.EQ.2.0) THEN
                  psin2 = psin1 - osmnode(ifit)%fit_psin(1)
                ELSE
                  psin2 = psin1
                ENDIF
c                WRITE(logfp,*) 'PSIN1,2:',psin1,psin2
                SELECTCASE (NINT(osmnode(ifit)%fit_type))
                  CASE ( 1)  ! Polynomial
                    val = osmnode(ifit)%fit_p(1)            +
     .                    osmnode(ifit)%fit_p(2) * psin2    +
     .                    osmnode(ifit)%fit_p(3) * psin2**2 +
     .                    osmnode(ifit)%fit_p(4) * psin2**3 +
     .                    osmnode(ifit)%fit_p(5) * psin2**4 +
     .                    osmnode(ifit)%fit_p(6) * psin2**5
                  CASE ( 2)  ! Exponential
                    val = osmnode(ifit)%fit_p(1) * 
     .                    EXP(osmnode(ifit)%fit_p(2)*psin2) + 
     .                    osmnode(ifit)%fit_p(3)
                  CASE ( 3)  ! TANH 
                    p(0:4) = osmnode(ifit)%fit_p(1:5)
                    v0 = (p(0) - psin2) / (2.0 * p(1))               ! from /home/mastts/lib/edgefunctionats.pro
                    v1 = ((1. + p(3) * v0) * EXP(v0) - EXP(-v0)) /   ! from /home/mastts/bck/mtanh.pro (right one?) 
     .                                      (EXP(v0) + EXP(-v0))
                    val = (p(2) - p(4)) / 2.0 * (v1 + 1.0) + p(4)
                  CASEDEFAULT
                    CALL ER('AssignNodeValues_2',
     .                      'Unknown data type',*99)
                ENDSELECT
                SELECTCASE (NINT(osmnode(ifit)%fit_quantity))
                  CASE (1)  
                    IF (osmnode(ifit)%fit_type.EQ.3.0)  ! Special for TANH fit from ts.pro
     .                val = val * 1.0E+19  
                    ne(inode) = val
                  CASE (2)  
                    vb(inode) = val 
                  CASE (4)  
                    te(inode) = val 
                  CASEDEFAULT
                    CALL ER('AssignNodeValues','Unknown data type',*99)
                ENDSELECT
              ENDIF
            ENDDO
c         --------------------------------------------------------------
          CASE (6)
c...        Core and pedestal fits based on shape and engineering
c           parameters:
            nfit = 0
            DO i4 = i1+1, osmnnode
              IF (osmnode(i4)%type.NE.0.0) EXIT  
              nfit = nfit + 1
            ENDDO
            WRITE(logfp,*) 'NFIT 6:',nfit
            DO i4 = 1, nfit
              ifit = i1 + i4
              WRITE(logfp,*) 'FIT:',val
              WRITE(logfp,*) '   :',osmnode(ifit)%fit_type
              WRITE(logfp,*) '   :',osmnode(ifit)%fit_quantity
              WRITE(logfp,*) '   :',osmnode(ifit)%fit_p(1:8)
              IF     (osmnode(ifit)%fit_type.EQ. 0.0) THEN 
              ELSEIF (osmnode(ifit)%fit_type.EQ.-1.0) THEN 
                val = val - osmnode(ifit)%fit_quantity  ! Hack
              ELSE
                type = NINT(osmnode(ifit)%fit_type)
                osmnode(ifit)%tube_range = osmnode(i0)%tube_range
                SELECTCASE (type)
                  CASE (1:2)  ! Core + pedestal + SOL
                    CALL osm_UpstreamProfile(type,ifit,val,coord,result)
                  CASEDEFAULT
                    CALL ER('AssignNodeValues_2','Unknown fit '//
     .                      'type for MODE=6',*99)
                ENDSELECT
                SELECTCASE (NINT(osmnode(ifit)%fit_quantity))
                  CASE (1)  
                    ne(inode) = result(1)
                  CASE (2)  
                    STOP 'NO VB DATA YET'
                  CASE (4)  
                    te(inode) = result(2)
                    ti(inode) = result(3)
                  CASE (5)  
                    ti(inode) = result(3)
                  CASEDEFAULT
                    CALL ER('AssignNodeValues_2','Unknown data '//
     .                      'assignment for MODE=6',*99)
                ENDSELECT
              ENDIF
            ENDDO                
c         ----------------------------------------------------------
          CASE (7)
c...        Exponential decay for v,T (p not allowed), using v_perp and L for n:
            C = expon
            A = v0 - v1
            B = v1
            IF (vc) vb(inode) = A * EXP(-val / C) + B
            A = te0 - te1
            B = te1 
            IF (tec) te(inode) = A * EXP(-val / C) + B
            A = ti0 - ti1
            B = ti1 
            IF (tic) ti(inode) = A * EXP(-val / C) + B

            IF (pc)
     .        CALL ER('AssignNodeValues_New','Not allowed to '//
     .                'specify pressure with MODE=7',*99)
            IF (te(inode).EQ.0.0)
     .        CALL ER('AssignNodeValues_New','Te needs to be set '//
     .                'when calculating convective densities',*99)
            te_cs = te(inode)
            ti_cs = ti(inode)
            IF (ti(inode).EQ.0.0) ti_cs = te_cs * opt%ti_ratio(LO)
            IF (opt%ti_ratio(LO).NE.opt%ti_ratio(HI))
     .        CALL ER('AssignNodeValues_New','Ti/Te ratio poorly '//
     .                'defined, need to decide what to do',*99)

            SELECTCASE (opt%radvel)
              CASE (0)
                CALL ER('AssignNodeValues_New','Radial velocity '//
     .                  'option not set',*99) 
              CASE (1)
                radvel = opt%radvel_param(1)
              CASE DEFAULT
                CALL ER('AssignNodeValues_New','Invalid radial '//
     .                  'velocity option',*99) 
            ENDSELECT

            IF (ABS(expon).LT.2.0E-7) expon = 1.0
            C = tube(it)%smax * radvel / GetCs2(te_cs,ti_cs) ! * expon
c            C = tube(it)%smax * 10.0 / GetCs2(te_cs,ti_cs) ! * expon
c            C = tube(it)%smax * 30.0 / GetCs2(te_cs,ti_cs) ! * expon
c            C = tube(it)%smax * 100.0 / GetCs2(te_cs,ti_cs) ! * expon
            A = n0
            B = 0.0
            IF (nc) THEN
              ne(inode) = A * EXP(-val / C) + B
              IF (debug) THEN 
                WRITE(logfp,*) '  L    :',tube(it)%smax
                WRITE(logfp,*) '  Te   :',te_cs
                WRITE(logfp,*) '  Ti   :',ti_cs
                WRITE(logfp,*) '  cs   :',GetCs2(te_cs,ti_cs)
                WRITE(logfp,*) '  A,B,C:',A,B,C
                WRITE(logfp,*) '  val  :',val
                WRITE(logfp,*) '  expon:',expon
                WRITE(logfp,*) '  ne   :',ne(inode)
              ENDIF
              ne(inode) = MAX(ne(inode),1.0E+14)
            ENDIF
          CASE DEFAULT
            CALL ER('AssignNodeValues_2','Invalid MODE B',*99)   
        ENDSELECT

c...    Store node values:
        IF (node_valid) THEN 
          node_n = node_n + 1
          node_i(node_n) = index
          node_s(node_n)%s     = s (inode)
          node_s(node_n)%ne    = ne(inode)
          node_s(node_n)%v     = vb(inode)
          node_s(node_n)%pe    = pe(inode)
          node_s(node_n)%te    = te(inode)
          node_s(node_n)%ti(1) = ti(inode)            

          node_s(node_n)%par_mode = osmnode(i3)%par_mode
          node_s(node_n)%par_exp  = osmnode(i3)%par_exp
          node_s(node_n)%par_set  = osmnode(i3)%par_set
        ENDIF

        IF (logop.GT.0) THEN
          WRITE(logfp,*) 
          DO i4 = 1, node_n
            WRITE(logfp,'(A,3I6,E10.2,3E10.2,2F10.2)') 
     .        ' >>> BUILDING NODES:',i4,node_i(i4),
     .        node_s(i4)%icell,node_s(i4)%s,
     .        node_s(i4)%ne,
     .        node_s(i4)%v,
     .        node_s(i4)%pe,
     .        node_s(i4)%te,node_s(i4)%ti(1)
          ENDDO
        ENDIF

      ENDDO  ! END OF OSMNMONE LOOP
c     ------------------------------------------------------------------

c...  Check for target data:
      DO i1 = 1, node_n
        IF (node_s(i1)%par_mode.EQ.-1) THEN
          node_s(node_n+1) = node_s(1      )
          node_i(node_n+1) = node_i(1       )
          node_s(1       ) = node_s(i1      )
          node_i(1       ) = node_i(i1      )
          node_s(i1      ) = node_s(node_n+1)
          node_i(i1      ) = node_i(node_n+1)
          node_s(1)%par_set = 2          
          EXIT
        ENDIF
      ENDDO
      IF (i1.EQ.node_n+1) THEN
        node_n = node_n + 1
        DO i1 = node_n, 2, -1
          node_i(i1) = node_i(i1-1)
          node_s(i1) = node_s(i1-1)
        ENDDO
        node_i(1) = 0
        node_s(1)%ne = 0.0
        node_s(1)%v  = 0.0
        node_s(1)%pe = 0.0
        node_s(1)%te = 0.0
        node_s(1)%ti(1) = 0.0
        node_s(1)%par_mode = 0
        node_s(1)%par_exp  = 0.0
        node_s(1)%par_set  = 0
      ENDIF
      DO i1 = 1, node_n
        IF (node_s(i1)%par_mode.EQ.-2) THEN
          node_s(node_n+1) = node_s(i1      )
          node_i(node_n+1) = node_i(i1      )
          node_s(i1      ) = node_s(node_n  )
          node_i(i1      ) = node_i(node_n  )
          node_s(node_n  ) = node_s(node_n+1)
          node_i(node_n  ) = node_i(node_n+1)
          node_s(node_n)%par_set = 2          
          EXIT
        ENDIF
      ENDDO
      IF (i1.EQ.node_n+1) THEN
        node_n = node_n + 1
        node_i(node_n) = 0
        node_s(node_n)%ne = 0.0
        node_s(node_n)%v  = 0.0
        node_s(node_n)%pe = 0.0
        node_s(node_n)%te = 0.0
        node_s(node_n)%ti(1) = 0.0
        node_s(node_n)%par_mode = 0
        node_s(node_n)%par_exp  = 0.0
        node_s(node_n)%par_set  = 0
      ENDIF
c...  Define target node location defaults:
      node_s(1     )%s = 0.0
      node_s(node_n)%s = tube(it)%smax
c      node_s(1     )%icell = 0
c      node_s(node_n)%icell = tube(it)%n + 1
c...  Assign cell indices:
       WRITE(logfp,*) 'length',tube(it)%smax
      DO i1 = 1, node_n
c      DO i1 = 2, node_n-1
        IF     (node_s(i1)%s.LT.0.00001) THEN
          node_s(i1)%s     = 0.0
          node_s(i1)%icell = 0
        ELSEIF (node_s(i1)%s.GT.0.99999*tube(it)%smax) THEN
          node_s(i1)%s     = tube(it)%smax
          node_s(i1)%icell = tube(it)%n + 1
        ELSE
          DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)
            IF (node_s(i1)%s.GE.cell(ic)%sbnd(1).AND.
     .          node_s(i1)%s.LE.cell(ic)%sbnd(2)) THEN
              node_s(i1)%icell = ic - tube(it)%cell_index(LO) + 1
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      IF (logop.GT.0) THEN
        WRITE(logfp,*) 
        DO i1 = 1, node_n
          WRITE(logfp,'(A,3I6,F10.2,3E10.2,2F10.2,5X,I6)') 
     .      'NODE TARGETS:',i1,node_i(i1),
     .      node_s(i1)%icell,node_s(i1)%s,
     .      node_s(i1)%ne,
     .      node_s(i1)%v,
     .      node_s(i1)%pe,
     .      node_s(i1)%te,node_s(i1)%ti(1), tube(it)%n
        ENDDO
      ENDIF

c...  Check symmetry node is specified:
      tube2(itube)%state = ibclr(tube2(itube)%state,0)
      i2 = 0
      DO i1 = 1, node_n
        IF (node_i(i1).EQ.2) i2 = i2 + 1
      ENDDO
      IF     (i2.EQ.0) THEN
c...    None found so add a symmetry node proxy:
        node_n = node_n + 1
        node_i(node_n) = 2 
        node_s(node_n)%type     = 2.0
        node_s(node_n)%par_mode = 1
        node_s(node_n)%par_exp  = 1.0
        node_s(node_n)%par_set  = 2
        node_s(node_n)%s = 0.5 * tube(it)%smax
c...    Assign cell index:
        DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)
          IF (node_s(node_n)%s.GE.cell(ic)%sbnd(1).AND.
     .        node_s(node_n)%s.LE.cell(ic)%sbnd(2))THEN
            node_s(node_n)%icell = ic - tube(it)%cell_index(LO) + 1
          ENDIF
        ENDDO
        IF (it.LT.grid%isep) THEN
          node_s(node_n)%ne    = 1.0E+19
          node_s(node_n)%te    = 100.0
          node_s(node_n)%ti(1) = 200.0
        ELSE
          node_s(node_n)%ne    = 1.0E+18
          node_s(node_n)%te    = 10.0
          node_s(node_n)%ti(1) = 20.0
        ENDIF
        tube2(itube)%state = ibset(tube2(itube)%state,0)
        WRITE(logfp,*) 'default symmetry point applied'
        WRITE(dummy,'(A,I6,A)') 'Symmetry node not identified for '//
     .                          'ITUBE = ',it,', applying default'
        IF (ntube.GT.200) THEN  ! Trying to avoid clutter on the terminal window 
          WRITE(PINOUT,*) TRIM(dummy)
          IF (default_message) THEN
            WRITE(0,*) 'WARNING AssignNodes: Default values applied'
            default_message = .FALSE.
          ENDIF
        ELSEIF (.NOT.suppress_screen) THEN
          CALL WN('AssignNodeValues_New',TRIM(dummy))
        ENDIF
      ELSEIF (i2.GT.1) THEN
c...    Try to combine nodes if there is more than one:
        DO i1 = 1, node_n
          IF (node_i(i1).EQ.2) EXIT
        ENDDO          
c...    Check that all the symmetry nodes are at the same location:
        DO i2 = i1+1, node_n
          IF (node_i(i2).EQ.2.AND.
     .        node_s(i1)%icell.NE.node_s(i2)%icell) 
     .        CALL ER('AssignNodeValues_New','Multiple '//
     .                'symmetry points not at same location',*99)
        ENDDO
        i2 = i1
        DO WHILE(i2.LT.node_n)
          i2 = i2 + 1
          IF (node_i(i2).EQ.2) THEN
c...       Check for quantity assigned multiple values:
           IF ((node_s(i1)%ne.NE.0.0.AND.node_s(i2)%ne.NE.0.0).OR.
     .         (node_s(i1)%v .NE.0.0.AND.node_s(i2)%v .NE.0.0).OR.
     .         (node_s(i1)%pe.NE.0.0.AND.node_s(i2)%pe.NE.0.0).OR.
     .         (node_s(i1)%te.NE.0.0.AND.node_s(i2)%te.NE.0.0)) 
     .       CALL ER('AssignNodeValues_New','Degenerate '//
     .               'symmetry point node found',*99)
c...       Copy over data:
           IF (node_s(i2)%ne.NE.0.0) node_s(i1)%ne = node_s(i2)%ne
           IF (node_s(i2)%v .NE.0.0) node_s(i1)%v  = node_s(i2)%v 
           IF (node_s(i2)%pe.NE.0.0) node_s(i1)%pe = node_s(i2)%pe
           IF (node_s(i2)%te.NE.0.0) node_s(i1)%te = node_s(i2)%te
           IF (node_s(i2)%ti(1).NE.0.) node_s(i1)%ti(1)=node_s(i2)%ti(1)
c...       Allow partial assignment of sheath-limited tubes by combining blocks:
           IF (node_s(i1)%par_set.EQ.0.AND.node_s(i2)%par_set.NE.0)
     .       node_s(i1)%par_set = node_s(i2)%par_set
c...       Delete degenerate node:
           IF (i2.LT.node_n) THEN
             node_s(i2:node_n-1) = node_s(i2+1:node_n)
             node_i(i2:node_n-1) = node_i(i2+1:node_n)
           ENDIF
           node_n = node_n - 1
          ENDIF
        ENDDO
      ENDIF





      IF (logop.GT.0) THEN
        WRITE(logfp,*) 
        DO i1 = 1, node_n
          WRITE(logfp,'(A,3I6,F10.2,3E10.2,2F10.2)') 
     .      'NODE SYMMETRY:',i1,node_i(i1),
     .      node_s(i1)%icell,node_s(i1)%s,
     .      node_s(i1)%ne,
     .      node_s(i1)%v,
     .      node_s(i1)%pe,
     .      node_s(i1)%te,node_s(i1)%ti(1)
        ENDDO
      ENDIF


c     THE ABOVE CODE FOR SPECIFICALLY CHECKING SYMMETRY POINTS MAY BE
c     REDUNDANT -- NEED TO CHECK

c...  Try to combine nodes if there is more than one:
      DO i1 = 1, node_n
        i2 = i1
        DO WHILE(i2.LT.node_n)
          i2 = i2 + 1
          IF (node_s(i1)%icell.EQ.node_s(i2)%icell) THEN
c...       Check for quantity assigned multiple values:
           IF ((node_s(i1)%ne.NE.0.0.AND.node_s(i2)%ne.NE.0.0).OR.
     .         (node_s(i1)%v .NE.0.0.AND.node_s(i2)%v .NE.0.0).OR.
     .         (node_s(i1)%pe.NE.0.0.AND.node_s(i2)%pe.NE.0.0).OR.
     .         (node_s(i1)%te.NE.0.0.AND.node_s(i2)%te.NE.0.0).OR.
     .         (node_s(i1)%ti(1).NE.0.0.AND.node_s(i2)%ti(1).NE.0.0)) 
     .       CALL ER('AssignNodeValues_New','Degenerate '//
     .               'symmetry point node found',*99)
c...       Copy over data:
           IF (node_s(i2)%ne.NE.0.0) node_s(i1)%ne = node_s(i2)%ne
           IF (node_s(i2)%v .NE.0.0) node_s(i1)%v  = node_s(i2)%v 
           IF (node_s(i2)%pe.NE.0.0) node_s(i1)%pe = node_s(i2)%pe
           IF (node_s(i2)%te.NE.0.0) node_s(i1)%te = node_s(i2)%te
           IF (node_s(i2)%ti(1).NE.0.) node_s(i1)%ti(1)=node_s(i2)%ti(1)
c...       Allow partial assignment of sheath-limited tubes by combining blocks:
           IF (node_s(i1)%par_set.EQ.0.AND.node_s(i2)%par_set.NE.0)
     .       node_s(i1)%par_set = node_s(i2)%par_set
c...       Delete degenerate node:
           IF (i2.LT.node_n) THEN
             node_s(i2:node_n-1) = node_s(i2+1:node_n)
             node_i(i2:node_n-1) = node_i(i2+1:node_n)
           ENDIF
           node_n = node_n - 1
          ENDIF
        ENDDO
      ENDDO






c...  Sort nodes based on s-distance along the field line:
      DO i1 = 1, node_n-1
        DO i2 = i1+1, node_n
          IF (node_s(i2)%s.LT.node_s(i1)%s) THEN
            node_s(0)  = node_s(i1)
            node_i(0)  = node_i(i1)
            node_s(i1) = node_s(i2)
            node_i(i1) = node_i(i2)
            node_s(i2) = node_s(0)
            node_i(i2) = node_i(0)
          ENDIF
c...      Also check that there isn't more than one node in each cell:
          IF (node_s(i1)%icell.EQ.node_s(i2)%icell) 
     .      CALL ER('AssignNodeValues_New','More than one node '//
     .              'detected in cell',*99)
        ENDDO
      ENDDO

      IF (logop.GT.0) THEN
        DO i1 = 1, node_n
          WRITE(logfp,'(A,3I6,F10.2,3E10.2,2F10.2)') 
     .      'NODE:',i1,node_i(i1),
     .      node_s(i1)%icell,node_s(i1)%s,
     .      node_s(i1)%ne,node_s(i1)%v,node_s(i1)%pe,
     .      node_s(i1)%te,node_s(i1)%ti(1)
        ENDDO
      ENDIF

c...  Set node indeces:
      DO i1 = 1, node_n 
        IF (node_i(i1).EQ.2) mnode = i1
      ENDDO
      nnode = node_n

c...  Assign values to nodes:
      node(1:nnode) = node_s(1:nnode)

      IF (logop.GT.0) THEN
        WRITE(logfp,*) 
        WRITE(logfp,'(A,2I6)') 'NODE A -:',nnode,mnode
        DO i1 = 1, node_n
          WRITE(logfp,'(A,4I6,F10.2,1P,3E10.2,0P,F10.3,
     .                  1P,E10.2,0P,2F10.2)') 
     .      'NODE -:',i1,node_i(i1),node(i1)%par_set,
     .      node(i1)%icell,node_s(i1)%s,
     .      node(i1)%jsat(1),node(i1)%ne,
     .      node(i1)%v,
     .      node(i1)%machno,
     .      node(i1)%pe,
     .      node(i1)%te,node(i1)%ti(1)
        ENDDO
      ENDIF

c...  Assign other quantites:
      node(1:nnode)%jsat(1) = 0.0
      node(1:nnode)%ni(1)   = 0.0
      node(1:nnode)%pi(1)   = 0.0
      DO i1 = 1, nnode
        IF (node(i1)%ti(1).EQ.0.0.AND.i1.LT.mnode) 
     .    node(i1)%ti(1) = node(i1)%te * opt%ti_ratio(LO) 
        IF (node(i1)%ti(1).EQ.0.0.AND.i1.EQ.mnode) 
     .    node(i1)%ti(1) = node(i1)%te * MAX(opt%ti_ratio(LO),
     .                                       opt%ti_ratio(HI))
        IF (node(i1)%ti(1).EQ.0.0.AND.i1.GT.mnode) 
     .    node(i1)%ti(1) = node(i1)%te * opt%ti_ratio(HI) 
      ENDDO
c      node(1      :mnode)%ti(1) =node(1      :mnode)%te*opt%ti_ratio(LO) 
c      node(mnode+1:nnode)%ti(1) =node(mnode+1:nnode)%te*opt%ti_ratio(HI) 
      node(1:nnode)%machno    = 0.0
      node(1:nnode)%epot      = 0.0
      node(1:nnode)%efield    = 0.0

c      IF (node(1    )%ne   .EQ.0.) node(1    )%ne   =tube(it)%ne(LO)
c      IF (node(nnode)%ne   .EQ.0.) node(nnode)%ne   =tube(it)%ne(HI)
c      IF (node(1    )%ti(1).EQ.0.) node(1    )%ti(1)=tube(it)%ti(LO,1)
c      IF (node(nnode)%ti(1).EQ.0.) node(nnode)%ti(1)=tube(it)%ti(HI,1)


c...  Set electron pressure from a reference node:
      DO i1 = mnode-1, 1, -1
        IF (node(i1)%pe.LT.-88.1.OR.node(i1)%pe.GT.-87.9) CYCLE
        node_pe = -1.0
        DO i2 = i1+1, mnode
          IF     (node(i2)%pe.GT.0.0) THEN
            node_pe = node(i2)%pe
          ELSEIF (node(i2)%ne.GT.0.0.AND.node(i2)%te.GT.0.0) THEN
            node_pe = node(i2)%ne * node(i2)%te
          ENDIF
        ENDDO
        IF (node_pe.EQ.-1.0) 
     .    CALL ER('AssignNodeValues_New','Pressure reference '//
     .            'requested but none found (LO)',*99) 
        node(i1)%pe = node_pe
      ENDDO
      DO i1 = mnode+1, nnode
        IF (node(i1)%pe.LT.-88.1.OR.node(i1)%pe.GT.-87.9) CYCLE
        node_pe = -1.0
        DO i2 = i1-1, mnode, -1
          IF     (node(i2)%pe.GT.0.0) THEN
            node_pe = node(i2)%pe
          ELSEIF (node(i2)%ne.GT.0.0.AND.node(i2)%te.GT.0.0) THEN
            node_pe = node(i2)%ne * (node(i2)%te + node(i2)%ti(ion))  ! changed 28.04.12
c            node_pe = node(i2)%ne * node(i2)%te
          ENDIF
        ENDDO
        IF (node_pe.EQ.-1.0) 
     .    CALL ER('AssignNodeValues_New','Pressure reference '//
     .            'requested but none found (HI)',*99) 
        node(i1)%pe = node_pe
      ENDDO


c...  Set electron pressure from a reference node:
      DO i1 = 1, nnode
       IF ((node(i1)%ne     .GT.-77.1.AND.node(i1)%ne     .LT.-76.9).OR.
     .     (node(i1)%v      .GT.-77.1.AND.node(i1)%v      .LT.-76.9).OR.
     .     (node(i1)%pe     .GT.-77.1.AND.node(i1)%pe     .LT.-76.9).OR.
     .     (node(i1)%te     .GT.-77.1.AND.node(i1)%te     .LT.-76.9).OR.
     .     (node(i1)%ti(ion).GT.-77.1.AND.node(i1)%ti(ion).LT.-76.9))
     .    THEN

          IF (ref_ntube.EQ.0.OR.ref_ntube.EQ.1) 
     .     CALL ER('AssignNodeValues_New','-77.0 found but reference'//
     .             'solution not available',*99)     

          ic1 = tube(itube)%cell_index(LO)
          ic2 = tube(itube)%cell_index(HI)     
          ALLOCATE(fluid_tmp(ic2-ic1+1,nion))
          tube_tmp = tube(itube)
          CALL InterpolateReferencePlasma(tube_tmp,nion,
     .           fluid_tmp(1:ic2-ic1+1,1:nion),cell(ic1:ic2))         
          EXIT
        ENDIF
      ENDDO
      IF (ALLOCATED(fluid_tmp)) THEN 
        DO i1 = 1, nnode
          IF (node(i1)%ne     .GT.-77.1.AND.node(i1)%ne     .LT.-76.9)
     .     node(i1)%ne       = fluid_tmp(node(i1)%icell,ion)%ne
          IF (node(i1)%v      .GT.-77.1.AND.node(i1)%v      .LT.-76.9) 
     .      THEN
            cs = GetCs2(fluid_tmp(node(i1)%icell,ion)%te,
     .                  fluid_tmp(node(i1)%icell,ion)%ti)
            node(i1)%machno  = fluid_tmp(node(i1)%icell,ion)%vi / cs 
            node(i1)%machno  = MIN(0.5,MAX(-0.5,node(i1)%machno))   ! *** TEMPORARY UNTIL THE SOLVER IS STABILIZED FOR HIGHER M ***
            node(i1)%v       = 0.0D0
          ENDIF
          IF (node(i1)%pe     .GT.-77.1.AND.node(i1)%pe     .LT.-76.9) 
     .      node(i1)%pe      = fluid_tmp(node(i1)%icell,ion)%ne *
     .                         fluid_tmp(node(i1)%icell,ion)%te  
          IF (node(i1)%te     .GT.-77.1.AND.node(i1)%te     .LT.-76.9)
     .      node(i1)%te      = fluid_tmp(node(i1)%icell,ion)%te
          IF (node(i1)%ti(ion).GT.-77.1.AND.node(i1)%ti(ion).LT.-76.9)
     .      node(i1)%ti(ion) = fluid_tmp(node(i1)%icell,ion)%ti
        ENDDO
        DEALLOCATE(fluid_tmp)
      ENDIF



c...  Setup target nodes:
      DO itarget = LO, HI
c...    Identify the target node and the nearest upstream node where dependence
c       of the target values on the upstream node has been identified:
        IF (itarget.EQ.LO) THEN
          i1 = 1
          DO i2 = 2, mnode
            IF (node(i2)%par_set.NE.0) EXIT
          ENDDO
          IF (i2.EQ.mnode+1) i2 = 2
        ELSE
          i1 = nnode
          DO i2 = nnode-1, mnode, -1
            IF (node(i2)%par_set.NE.0) EXIT
          ENDDO
          IF (i2.EQ.mnode-1) i2 = nnode - 1
        ENDIF
c...    Assign target values in TUBE(IT) based on the corresponding target node:
        SELECTCASE(node(i1)%par_set)
c         --------------------------------------------------------------
          CASE (0) 
c         --------------------------------------------------------------
          CASE (2) 
            IF (node(i1)%ne.NE.0.0) THEN
              IF (node(i1)%ne.LT.1.0E+10) THEN
                tube(it)%jsat(itarget,ion) = node(i1)%ne
              ELSE
                IF (node(i1)%te.NE.0.0) THEN
                  tube(it)%jsat(itarget,ion) = 
     .              GetJsat2(node(i1)%te,
     .                       node(i1)%ti(ion),
     .                       node(i1)%ne,1.0) 
                ELSE
                  tube(it)%jsat(itarget,ion) = 
     .              GetJsat2(tube(it)%te(itarget),
     .                       tube(it)%ti(itarget,ion),
     .                       node(i1)%ne,1.0) 
                ENDIF
              ENDIF
            ENDIF
            IF (node(i1)%te     .NE.0.0) 
     .        tube(it)%te(itarget    ) = node(i1)%te
            IF (node(i1)%ti(ion).NE.0.0) 
     .        tube(it)%ti(itarget,ion) = node(i1)%ti(ion)
c         --------------------------------------------------------------
          CASE DEFAULT
            CALL ER('_New','Unknown PAR_SET for target node',*99) 
c         --------------------------------------------------------------
        ENDSELECT
c...    Assign the target values in TUBE(IT) based on an upstream node:
        SELECTCASE(node(i2)%par_set)
c         --------------------------------------------------------------
          CASE (0) 
c         --------------------------------------------------------------
          CASE (1) 
            tube(it)%te(itarget)     = node(i2)%te
            tube(it)%ti(itarget,ion) = node(i2)%ti(ion)
c         --------------------------------------------------------------
          CASE (2:4) 

            IF (node(i2)%par_set.EQ.2.AND.node(i1)%te.NE.0.0) THEN 
              IF (te_warning) THEN
                WRITE(0,*) 'WARNING: PAR_SET=2 but target '//
     .           'Te data already assigned, overwriting'
                te_warning = .FALSE.
              ENDIF
              node(i1)%te = 0.0
            ENDIF

            IF (node(i1)%te.EQ.0.0) THEN
              SELECTCASE (node(i2)%par_set)
                CASE (2:3) 
                  node_te = node(i2)%te 
                CASE (4)  
                  IF (node_i(i2).NE.2)
     .              CALL ER('_New','Symmetry node required',*99) 
	    
                  IF (itarget.EQ.LO) THEN 
                    L = node(i2)%s
                  ELSE	          
                    L = tube(it)%smax - node(i2)%s
                  ENDIF
	    
                  nustar = (1.0E-16*node(i2)%ne) * L / node(i2)%te**2     
	    
                  frac = (nustar - 1.0) / (50.0 - 1.0)        ! parameter  *** THE  1.0 DOESN'T DO ANYTHING HERE! ***
c                  frac = (nustar - 10.0) / (50.0 - 10.0)        ! parameter
                  write(logfp,*) 'frac    =',frac
                  frac = MIN(MAX(0.0,frac),1.0)
	          
                  node_te = node(i2)%te
                  node_te = (1.0-frac) *     node_te + 
     .                           frac  * MIN(node_te,5.0)                  ! parameter
	    
                  IF (debug) THEN
                    write(logfp,*) 'tube  = ',it
                    write(logfp,*) 'node  = ',i2,node_i(i2)	          
                    write(logfp,*) 'sdat  = ',node(i2)%s,tube(it)%smax
                    write(logfp,*) 'n     = ',node(i2)%ne
                    write(logfp,*) 'te    = ',node(i2)%te
                    write(logfp,*) 'L     = ',L
                    write(logfp,*) 'nustar  =',nustar
                    write(logfp,*) 'frac    =',frac
                    write(logfp,*) 'node_te =',node_te,'  <---'
                  ENDIF
	    
              ENDSELECT
              tube(it)%te(itarget) = node_te
            ENDIF
	     
c            IF ((node(i2)%par_set.EQ.2                            ).OR.
c     .          (node(i2)%par_set.EQ.3.AND.node(i1)%ti(ion).EQ.0.0))
c     .        tube(it)%ti(itarget,ion) = node(i2)%ti(ion)


            IF (node(i2)%par_set.EQ.2.AND.node(i1)%ti(ion).NE.0.0) THEN
              IF (ti_warning) THEN
                WRITE(0,*) 'WARNING: PAR_SET=2 but target '//
     .           'Ti data already assigned, overwriting'
                ti_warning = .FALSE.
              ENDIF
              node(i1)%ti(ion) = 0.0
            ENDIF

            IF (node(i1)%ti(ion).EQ.0.0) THEN
              SELECTCASE (node(i2)%par_set)
                CASE (2:3) 
                  node_ti = node(i2)%ti(ion)
                CASE (4)  
                  IF (node_i(i2).NE.2)
     .              CALL ER('_New','Symmetry node required',*99) 
	    
                  IF (itarget.EQ.LO) THEN 
                    L = node(i2)%s
                  ELSE	          
                    L = tube(it)%smax - node(i2)%s
                  ENDIF
	    
                  nustar = (1.0E-16 * node(i2)%ne) * L / 
     .                     node(i2)%ti(ion)**2     
	          
                  frac = (nustar - 0.1) / (50.0 - 0.1)        ! parameter  *** THE  0.1 DOESN'T DO ANYTHING HERE! ***
                  frac = MIN(MAX(0.0,frac),1.0)
	          
                  node_ti = node(i2)%ti(ion)
                  node_ti = (1.0-frac) *     node_ti + 
     .                           frac  * MIN(node_ti,5.0)                  ! parameter
	    
                  IF (debug) THEN
                    write(logfp,*) 'ti      = ',node(i2)%ti(ion)
                    write(logfp,*) 'nustar  =',nustar
                    write(logfp,*) 'frac    =',frac
                    write(logfp,*) 'node_ti =',node_ti,'  <---'
                  ENDIF
	    
              ENDSELECT
              tube(it)%ti(itarget,ion) = node_ti
            ENDIF


            IF (node(i2)%par_set.EQ.2.AND.
     .          (node(i1)%ne.NE.0.0.OR.node(i1)%pe.NE.0.0)) THEN
              IF (ne_warning) THEN 
                WRITE(0,*) 'WARNING AssignNodeValues_New: ne and/or '//
     .                     'pe already assigned, overwriting'
                ne_warning = .FALSE.
              ENDIF
              node(i1)%ne = 0.0
              node(i1)%pe = 0.0
            ENDIF

            IF (node(i1)%ne.EQ.0.0.AND.node(i1)%pe.EQ.0.0)
     .        tube(it)%jsat(itarget,ion) = REAL(-i2)

c         --------------------------------------------------------------
          CASE DEFAULT
            CALL ER('_New','Unknown PAR_SET',*99) 
c         --------------------------------------------------------------
        ENDSELECT
c...    Assign target node values:
        node(i1)%jsat(ion) = tube(it)%jsat(itarget,ion)
        node(i1)%ne        = 0.0
        node(i1)%pe        = 0.0
        node(i1)%te        = tube(it)%te(itarget)
        node(i1)%ti(ion)   = tube(it)%ti(itarget,ion)
      ENDDO

c...  Sort velocities from Mach numbers:
      DO i1 = 1, node_n
        IF(node(i1)%v.NE.0.0.AND.ABS(node(i1)%v).LT.10.0) THEN
          node(i1)%machno = node(i1)%v
          node(i1)%v      = 0.0
        ENDIF
      ENDDO


      IF (logop.GT.0) THEN
        WRITE(logfp,*) 
        WRITE(logfp,'(A,2I6)') 'NODE C -:',nnode,mnode
        DO i1 = 1, node_n
          WRITE(logfp,'(A,4I6,F10.2,1P,3E10.2,0P,F10.3,
     .                  1P,E10.2,0P,2F10.2)') 
     .      'NODE -:',i1,node_i(i1),node(i1)%par_set,
     .      node(i1)%icell,node_s(i1)%s,
     .      node(i1)%jsat(1),node(i1)%ne,
     .      node(i1)%v,
     .      node(i1)%machno,
     .      node(i1)%pe,
     .      node(i1)%te,node(i1)%ti(1)

        ENDDO
      ENDIF

c...  Upscale the electron pressure in node%pe to total static pressure:
      DO i1 = 1, node_n 
        IF (node(i1)%pe.NE.0.0) THEN

          IF (node(i1)%te     .EQ.0.0) 
     .      CALL ER('AssignNodeValues_New','Attempting to re-scale '//
     .              'pressure but Te not defined',*99)
          IF (node(i1)%ti(ion).EQ.0.0) 
     .      CALL ER('AssignNodeValues_New','Attempting to re-scale '//
     .              'pressure but Ti not defined',*99)

          node(i1)%pe = node(i1)%pe * (node(i1)%te + node(i1)%ti(1)) /    ! changed 28.04.12
     .                                 node(i1)%te
        ENDIF
      ENDDO


      IF (logop.GT.0) THEN
        WRITE(logfp,*) 
        WRITE(logfp,'(A,2I6)') 'NODE B -:',nnode,mnode
        DO i1 = 1, node_n
          WRITE(logfp,'(A,4I6,F10.2,1P,3E10.2,0P,F10.3,
     .                  1P,E10.2,0P,2F10.2)') 
     .      'NODE -:',i1,node_i(i1),node(i1)%par_set,
     .      node(i1)%icell,node_s(i1)%s,
     .      node(i1)%jsat(1),node(i1)%ne,
     .      node(i1)%v,
     .      node(i1)%machno,
     .      node(i1)%pe,
     .      node(i1)%te,node(i1)%ti(1)

        ENDDO

      ENDIF


c      IF (ALLOCATED(fluid_tmp)) STOP 'dfsdfsd'


      RETURN
99    WRITE(0,*) ' ITUBE=',itube
      WRITE(0,*) ' MNODE=',mnode
      WRITE(0,*) ' NNODE=',nnode
      WRITE(0,*) ' TYPE =',tube(itube)%type
      WRITE(0,*) ' MODE =',mode
      WRITE(0,*) ' COORD=',coord
      IF (i1.GE.1.AND.i1.LT.node_n) THEN
        WRITE(0,*) ' 1 NODE NE =',node_s(i1)%ne
        WRITE(0,*) ' 1 NODE V  =',node_s(i1)%v
        WRITE(0,*) ' 1 NODE PE =',node_s(i1)%pe
        WRITE(0,*) ' 1 NODE TE =',node_s(i1)%te
        WRITE(0,*) ' 1 NODE TI =',node_s(i1)%ti
      ENDIF
      IF (i2.GE.1.AND.i2.LT.node_n) THEN
        WRITE(0,*) ' 2 NODE NE =',node_s(i2)%ne
        WRITE(0,*) ' 2 NODE V  =',node_s(i2)%v
        WRITE(0,*) ' 2 NODE PE =',node_s(i2)%pe
        WRITE(0,*) ' 2 NODE TE =',node_s(i2)%te
        WRITE(0,*) ' 2 NODE TI =',node_s(i2)%ti
      ENDIF
      CALL DumpData_OSM('output.error_node','Trouble building nodes')      
      STOP
      END


