c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE SetTargetConditions
      USE mod_sol28_global
      IMPLICIT none


      CALL SetTargetConditions_Legacy


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
      SUBROUTINE UnzipFile(fname)
      IMPLICIT none

      CHARACTER fname*(*)

      INTEGER   status
      CHARACTER command*1024

      IF     (fname(LEN_TRIM(fname)-2:LEN_TRIM(fname)).EQ.'zip') THEN
        command = 'unzip '//TRIM(fname)
        fname(LEN_TRIM(fname)-3:LEN_TRIM(fname)) = ' '
      ELSEIF (fname(LEN_TRIM(fname)-1:LEN_TRIM(fname)).EQ.'gz' ) THEN
        command = 'gunzip '//TRIM(fname)
        fname(LEN_TRIM(fname)-2:LEN_TRIM(fname)) = ' '
      ENDIF

      CALL CIssue(TRIM(command),status)        
      IF (status.NE.0) CALL ER('UnzipFiles','Dismal failure',*99)

      RETURN
 99   WRITE(0,*) '  FILE NAME = ',TRIM(fname)
      WRITE(0,*) '  COMMAND   = ',TRIM(command)
      WRITE(0,*) '  ERROR     = ',status
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
      SUBROUTINE SetTargetConditions_Legacy
      USE mod_legacy
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

c.... OSM:
      INTEGER itube,ishift


c...  Local:
      REAL GetRelaxationFraction

      INTEGER iitersol

      INTEGER ir,ii,ind1,ind2,i1,i2,i3,region,method,ring,target

      LOGICAL specific,apply,repeat,cycle_loop
      REAL    dum1,dum2,dum3,dum4,exp2,exp3,exp4,dpsin,psin0,frac,
     .        tedat(MAXNRS),tidat(MAXNRS),nedat(MAXNRS)


      CALL AssignLegacyVariables

      IF (log.GT.0) WRITE(logfp,*) 'SETTING TARGET CONDITIONS'
      
      IF (tarninter(LO).EQ.0.OR.tarninter(HI).EQ.0) THEN
        CALL WN('InterpolateTargetData','No data to '//
     .          'interpolate?')
        RETURN
      ENDIF

      IF (connected) THEN
        WRITE(0,*) 'CHEATING ON PSITARG'
        DO ir = irsep, nrs
          psitarg(ir,1) = REAL(ir-irsep)/100.0 + 1.0
          psitarg(ir,2) = REAL(ir-irsep)/100.0 + 1.0
        ENDDO
      ENDIF

      WRITE(SLOUT,*)
      WRITE(SLOUT,*) 'INTERPOLATING TARGET DATA'

c...  Need to check that PSITARG is assigned:
      IF (psitarg(irsep,2).EQ.0.0)
     .  CALL ER('InterpolateTargetData','PSITARG does not appear to '//
     .                                  'be assigned',*99)

      repeat = .TRUE.


      frac = GetRelaxationFraction()

c...  Overrides any explicit assignments in the input file.  Not strictly 
c     required, but leaving this in until it is a problem:
      nlpdato = 0
      nlpdati = 0

c...  Assign target data (LPDATO, LPDATI) from the TARINTER
c     data listed in the DIVIMP/OEDGE input file:
 10   DO region = IKLO, IKHI
        IF (tarninter(region).EQ.0) CYCLE
        ii = 0
        DO ir = irsep, nrs

          ii = ii + 1

          IF (idring(ir).EQ.BOUNDARY) THEN
c...        Assign default values to the virtual rings:
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
            CYCLE
          ENDIF

c...      Check if target data has been assigned specifically for this ring: 
          specific = .FALSE.
          method = 1
          ind1 = 0
          ind2 = 0

          DO i1 = 1, tarninter(region)
            apply = .FALSE.
            IF (tarinter(i1,1,region).EQ.-1.0) THEN
              specific = .TRUE.

              IF     (ir.GE.NINT(tarinter(i1,2,region)).AND.
     .                ir.LE.NINT(tarinter(i1,3,region))) THEN
c...            Apply target data to ring:
                apply = .TRUE.
              ELSEIF (tarinter(i1,2,region).LT.0.0.AND.
     .                tarinter(i1,3,region).LT.0.0) THEN
c...            Make sure that ring group exists:
                IF (-NINT(tarinter(i1,3,region)).GT.grdntreg(region))
     .            CALL ER('InterpolateTargetData','Specified '//
     .                    'region index is not valid',*99)
c...            Check if ring IR is a member of the specified group(s)
c               of rings:
                DO i2 = -NINT(tarinter(i1,2,region)), 
     .                  -NINT(tarinter(i1,3,region))
                  DO i3 = 1, grdntseg(i2,region)
                    IF (ir.EQ.grdtseg(i3,i2,region)) apply = .TRUE. 
                  ENDDO
                ENDDO                           
              ENDIF
              IF (apply) THEN
c...            Interpolation scheme:
                IF (tarinter(i1,4,region).EQ.1.0) THEN              
c...              Versus PSIn, exponential decay:
                  method = 2
                ELSE
c...              Standard:
                  method = 1
                ENDIF
                ind1 = 0
                ind2 = tarninter(region)
                DO i2 = i1+1, tarninter(region)-1
                  IF (ind1.EQ.0.AND.
     .                tarinter(i2  ,1,region).NE.-1.0) ind1 = i2
                  IF (ind1.NE.0.AND.
     .                tarinter(i2+1,1,region).EQ.-1.0) THEN
                    ind2 = i2 
                    EXIT
                  ENDIF
                ENDDO
                WRITE(PINOUT,*) 'APPLYING:',ir,ind1,ind2,region
                IF (ind2-ind1+1.LT.2) 
     .            CALL ER('InterpolateTargetData','Insufficient '//
     .                    'target data for interpolation',*99)
              ENDIF

            ENDIF
          ENDDO

          IF (.NOT.specific) THEN  
c...        Target data was not specified for specific rings, so treat
c           the entire TARINTER array as target data to be interpolated:
            ind1 = 1
            ind2 = tarninter(region)
          ELSEIF (ind1.EQ.0.OR.ind2.EQ.0) THEN
c...        Data was not found for this ring:
            WRITE(0,*) 'WARNING: TARGET DATA NOT FOUND FOR',ir,region
            CYCLE
          ENDIF

          IF     (method.EQ.1) THEN
c...        Linearly interpolate target data from TARINTER arrays:
            IF (region.EQ.IKLO) THEN
              dum1 = psitarg(ir,2)
             ELSE
              dum1 = psitarg(ir,1)
            ENDIF

c...        Make a list:            
            i3 = 0
            DO i2 = ind1, ind2
              i3 = i3 + 1
              tedat(i3) = (1.0 - frac) * tarinter(i2,2,region) + 
     .                           frac  * tarinter(i2,5,region)
              tidat(i3) = (1.0 - frac) * tarinter(i2,3,region) + 
     .                           frac  * tarinter(i2,6,region)
              nedat(i3) = (1.0 - frac) * tarinter(i2,4,region) + 
     .                           frac  * tarinter(i2,7,region)
            ENDDO
            CALL Fitter(ind2-ind1+1,tarinter(ind1,1,region),tedat,
     .                  1,dum1,dum2,'LINEAR')
            CALL Fitter(ind2-ind1+1,tarinter(ind1,1,region),tidat,
     .                  1,dum1,dum3,'LINEAR')
            CALL Fitter(ind2-ind1+1,tarinter(ind1,1,region),nedat,
     .                  1,dum1,dum4,'LINEAR')
          ELSEIF (method.EQ.2) THEN
c...        Exponential decay along the target versus PSIn:
            target = 0
            IF (tarinter(ind1,1,region).LT.0.0) THEN
c...          Decay is scaled by the values for a specified target segement:
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
c...        Interpolate target profile parameters (in case target boundary
c           relaxation is being used):
            DO i2 = 1, 3
              tedat(i2) = (1.0 - frac) * tarinter(ind1+i2-1,2,region) +
     .                           frac  * tarinter(ind1+i2-1,5,region)
              tidat(i2) = (1.0 - frac) * tarinter(ind1+i2-1,3,region) +
     .                           frac  * tarinter(ind1+i2-1,6,region)
              nedat(i2) = (1.0 - frac) * tarinter(ind1+i2-1,4,region) +
     .                           frac  * tarinter(ind1+i2-1,7,region)
            ENDDO
c...        Allow over-ride (RELMODE=8) of profile parameters at each step,
c           when boundary condition relaxation is being used, REL_OPT=2,3: 
            IF (tedat(1).EQ.-1.0) THEN
              IF (relmode.EQ.8.AND.relexpdat(region,1,ir).NE.0.0) THEN
                tedat(1) = relexpdat(region,1,ir)
              ELSEIF (rel_step.EQ.0) THEN
c...            The 0th EIRENE iteration is (currently) being avoided when target relaxation
c               is being used, so just assign some arbitrary nedat(1)!=-1.0.  The value from 
c               RELEXPDAT will be assigned on a subsequent call to this routine before the 
c               plasma is calculated (when REL_STEP=1, see in BgPlasma routine):
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

c...        Set exponents:
            exp2 = EXP(-dpsin / tedat(3))
            exp3 = EXP(-dpsin / tidat(3))
            exp4 = EXP(-dpsin / nedat(3))
c...        Calculate target boundary conditions:
            dum2 = (tedat(1) - tedat(2)) * exp2 + tedat(2)
            dum3 = (tidat(1) - tidat(2)) * exp3 + tidat(2)
            dum4 = (nedat(1) - nedat(2)) * exp4 + nedat(2)
c...        (Need to replace the target data specifier which was overwritten 
c            above, nasty):
            IF (target.NE.0) tarinter(ind1,2,region) = REAL(target)
          ENDIF
c...      Assign LPDATx arrays (target data used to assign target data
c         KxDS arrays):
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
c...      End of IR loop:
        ENDDO
c...    End of REGION loop:
      ENDDO

      IF (repeat) THEN
c...    Need to do this twice:
        repeat = .FALSE.
        GOTO 10
      ENDIF

c      CALL Outputdata(85,'sdfds')
c      STOP 'sdgsdg'
	





c...  Copy data to OSM arrays:
      DO itube = 1, ntube
        IF (itube.LT.grid%isep) CYCLE

        cycle_loop = .TRUE.
        DO i1 = 1, opt%sol_n
          IF (itube.GE.opt%sol_tube(1,i1).AND.
     .        itube.LE.opt%sol_tube(2,i1)) cycle_loop = .FALSE.
        ENDDO
        IF (cycle_loop) CYCLE

c        WRITE(0,*) 'WHAT?',itube,opt%sol_n
c        WRITE(0,*) 'WHAT?',i1,opt%sol_tube(1:2,1)
c        WRITE(0,*) 'WHAT?',i1,opt_iteration(1)%sol_tube(1:2,1)

        DO ii = 1, nlpdato
c          ishift = 1
c          IF (itube.GE.grid%ipfz) ishift = 2
c           WRITE(0,*) 'LPDATI:',itube,tube(itube)%ir
          IF (NINT(lpdato(ii,1)).EQ.tube(itube)%ir) THEN
c            WRITE(0,*) '      :',ii,lpdato(ii,2)
            tube(itube)%te  (LO)   = lpdato(ii,2)
            tube(itube)%ti  (LO,1) = lpdato(ii,3)
            tube(itube)%jsat(LO,1) = lpdato(ii,4)
            WRITE(logfp,*) 'JSAT8:',itube,lpdato(ii,4)
          ENDIF
        ENDDO

        DO ii = 1, nlpdati
c          ishift = 1
c          IF (itube.GE.grid%ipfz) ishift = 2
          IF (NINT(lpdati(ii,1)).EQ.tube(itube)%ir) THEN
c          IF (NINT(lpdati(ii,1)).EQ.itube+ishift) THEN
            tube(itube)%te  (HI)   = lpdati(ii,2)
            tube(itube)%ti  (HI,1) = lpdati(ii,3)
            tube(itube)%jsat(HI,1) = lpdati(ii,4)
            WRITE(logfp,*) 'JSAT7:',itube,lpdati(ii,4)
          ENDIF
        ENDDO



      ENDDO


      RETURN
 99   WRITE(0,*) '  I!,REGION       =',i1,region
      WRITE(0,*) ' -NINT(TARINTER)  =',-NINT(tarinter(i1,3,region))
      WRITE(0,*) '  GRDNTREG(REGION)=',grdntreg(region)
      STOP
      END



c
c ======================================================================
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

      DO it = osmnode(ind1)%tube_range(1), 
     .        osmnode(ind1)%tube_range(2)
        IF (tube(it)%type.EQ.GRD_BOUNDARY) CYCLE

        DO i1 = ind0+1, ind1

          a1 = DBLE(osmnode(i1-1)%rad_x)
          a2 = DBLE(osmnode(i1-1)%rad_y)
          b1 = DBLE(osmnode(i1  )%rad_x)
          b2 = DBLE(osmnode(i1  )%rad_y)
          dist = REAL(DSQRT((a1 - b1)**2 + (a2 - b2)**2))

          DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)

c...        Assumed 1:1 mapping between grid and data:
            iobj = ic
            isrf = ABS(obj(iobj)%iside(1))
            ivtx(1:2) = srf(isrf)%ivtx(1:2)
            c1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
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
      DO it = osmnode(ind1)%tube_range(1), 
     .        osmnode(ind1)%tube_range(2)
        IF (tube(it)%type.EQ.GRD_BOUNDARY.OR.ic1(it).EQ.0) CYCLE

c        WRITE(0,*) 'PICKENS:',ir,ik1(ir),dist1(ir)
        IF     (it.EQ.itgive) THEN
          iccell(1) = ic1(it)
          itcell(1) = it
        ENDIF
        IF (clcell(2).EQ.0.0.OR.dist1(it).LT.clcell(2)) THEN
          iccell(2) = ic1(it)
          itcell(2) = it
          clcell(2) = dist1(it)
        ENDIF
        IF (clcell(3).EQ.0.0.OR.dist1(it).GT.clcell(3)) THEN
          iccell(3) = ic1(it)
          itcell(3) = it
          clcell(3) = dist1(it)
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
      SUBROUTINE LoadUpstreamData(fname,itube,coord,shift,
     .                            xcolumn,ycolumn,yval)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      CHARACTER, INTENT(IN)  :: fname*(*)
      INTEGER  , INTENT(IN)  :: itube,coord
      REAL     , INTENT(IN)  :: shift,xcolumn,ycolumn
      REAL     , INTENT(OUT) :: yval

      LOGICAL GetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

      INTEGER    MAXTDAT     ,MAXCOLS
      PARAMETER (MAXTDAT=1000,MAXCOLS=10)

      INTEGER   ndata,fp,i,j,n,ncolumns,xcol,ycol
      REAL      vdata(MAXTDAT,MAXCOLS),xval
      CHARACTER buffer*1024

c      WRITE(0,*) 'UPSTREAM: COORD2=',coord
c      WRITE(0,*) 'UPSTREAM: FILE='//TRIM(fname)//'<'

c...  Access data file:
      fp = 99
      OPEN(UNIT=fp,FILE=TRIM(fname),ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)
      DO WHILE (GetLine(fp,buffer,WITH_TAG))
c        WRITE(0,*) 'BUFFER:',TRIM(buffer)

c...    Isolate tag string:
        DO i = 2, LEN_TRIM(buffer)
          IF (buffer(i:i).EQ.'}') EXIT
        ENDDO

        n = LEN_TRIM(buffer)

c        WRITE(0,*) 'BUFFER >'//buffer(2:i-1)//'<'

        SELECTCASE (buffer(2:i-1))
          CASE ('NUMBER OF COLUMNS')
            READ(buffer(i+1:n),*) ncolumns
          CASE ('DATA LIST')
            j = 0
c              WRITE(0,*) 'NCOLUMNS:',ncolumns
            DO WHILE(GetLine(fp,buffer,NO_TAG))
              IF (buffer(1:1).EQ.'{') EXIT
              j = j + 1
c              WRITE(0,*) 'BUFFER -:',TRIM(buffer)
              READ(buffer,*) vdata(j,1:ncolumns)
            ENDDO
            ndata = j
          CASE ('END')
            EXIT
          CASE DEFAULT
            CALL ER('LoadUpstreamData','Unknown tag',*99)
        ENDSELECT

      ENDDO
      CLOSE(fp)

c... 
      SELECTCASE (coord)
        CASE (3)
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
        CALL Fitter(ndata,vdata(1,xcol),vdata(1,ycol),
     .              1,xval,yval,'LINEAR')

c        WRITE(88,*) 'FIT:',ndata
c        WRITE(88,*) '   :',vdata(1:ndata,xcol)
c        WRITE(88,*) '   :',vdata(1:ndata,ycol)
c        WRITE(88,*) '   :',xval,yval

      ENDIF

c      WRITE(0,*) 'XCOL,YCOL:',xcol,ycol
c      WRITE(0,*) 'SHIFT:',shift
c      WRITE(0,*) 'ITUBE,XVAL,YVAL:',itube,xval,yval

      RETURN
 98   WRITE(0,*) 'ERROR LoadUpstreamData: Data file not found'
      WRITE(0,*) '  FILE = ',TRIM(fname)
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE AssignNodeValues_New(itube,nnode,mnode,node,opt_tube)
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_geometry
      USE mod_legacy
      IMPLICIT none

      INTEGER, INTENT(IN)  :: itube
      INTEGER, INTENT(OUT) :: nnode,mnode       
      TYPE(type_node), INTENT(OUT) :: node(*)
      TYPE(type_options_osm) :: opt_tube

      REAL GetJsat2,GetRelaxationFraction

      INTEGER ic,ic1,it,it1,i0,i1,i2,i3,i4,ifit,index,mode,
     .        iobj,isrf,ivtx(2),nfit,icell(3),ion,coord
      LOGICAL nc,pc,tc,density,tetarget
      REAL    te(0:6),ne(0:6),s(0:6),pe(0:6),
     .        frac,t0,t1,n0,n1,A,B,C,expon,
     .        psin0,psin1,psin2,p0,p1,
     .        prb1,tmp1,val,val0,val1,val2,p(0:5),v0,v1,v2
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd


      INTEGER node_n,node_i(0:50)
      TYPE(type_node) :: node_s(0:50)


      node_n = 1

      frac = GetRelaxationFraction()

c...  Better/cleaner to pass the tube to this routine, and not need to
c     use mod_sol28_locals...?
      it = itube

      IF (osmnnode.EQ.0) 
     .  CALL ER('AssignNodeValues_2','Profiles not found',*99)


      DO i1 = 2, osmnnode

        te = 0.0
        pe = 0.0
        ne = 0.0
c...    
        i0 = i1 - 1

c        WRITE(0,*) 'I0,1:',i0,i1,osmnode(i0)%type,osmnode(i1)%type

        IF (osmnode(i0)%type.NE.osmnode(i1)%type.OR.
     .      osmnode(i1)%type.EQ.0.0.OR.
     .      osmnode(i1)%type.EQ.0.0) CYCLE

c...    Do not apply data if IT is outside specified range:
        IF (REAL(it).LT.osmnode(i1)%tube_range(1).OR.
     .      REAL(it).GT.osmnode(i1)%tube_range(2)) CYCLE


c...    Need to keep track of iteration currently being assigned, so
c       that scans in upstream with step can be done serially, rather than
c       having everything on the same line...




c...    Check that rings from different grid regions are not in the same group
c       of rings:
        DO i2 = osmnode(i1)%tube_range(1), osmnode(i1)%tube_range(2)-1
          IF (tube(i2)%type.NE.tube(i2+1)%type) THEN
            IF (log.GT.0.AND.tube(i2)%type.NE.GRD_CORE) THEN
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
c   5 - something strange...
c
c   coord = 1 - linear on line segment
c         = 2 - linear on line segment, but from first to last ring intersection
c         = 3 - PSIn 
c         = 4 - RHO
c         = 5 - PSIn (raw)

        index = NINT(osmnode(i1)%type)
        mode  = osmnode(i1)%rad_mode
        coord = osmnode(i1)%rad_coord
        expon = osmnode(i1)%rad_exp

c...    Decide if specified upstream data is density or pressure:
        density = .TRUE.   ! *** I DON'T LIKE THIS SOLUTION ***
        IF (index.EQ.3.AND.
     .      (tube(it)%type.EQ.GRD_SOL.AND..FALSE..OR.
     .       tube(it)%type.EQ.GRD_PFZ.AND..FALSE.))
     .    density = .FALSE.

        a1 = DBLE(osmnode(i1-1)%rad_x)
        a2 = DBLE(osmnode(i1-1)%rad_y)
        b1 = DBLE(osmnode(i1  )%rad_x)
        b2 = DBLE(osmnode(i1  )%rad_y)

        DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)

          s  = 0.0
          ne = 0.0
          pe = 0.0
          te = 0.0

c...      Assumed 1:1 mapping between grid and data:
          iobj = ic
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
     .        tcd.GE.0.0D0.AND.tcd.LT.1.0D0) THEN
c...        Intersecting between the line segment and the ring is found:

c...        These are here in case psin0 gets picked up when linking exponential
c           decay data to neighbouring ring:
            IF (tube(it)%type.EQ.GRD_PFZ) THEN
              psin0 = -RHI
              psin1 =  RHI
            ELSE
              psin0 =  RHI
              psin1 =  RHI
            ENDIF

            IF (.FALSE.) THEN
c...          Need a check to know if the ring is over-constrained:
            ELSEIF (index.GE.1.AND.index.LE.4) THEN

              s(index) = cell(ic)%sbnd(1) + SNGL(tcd) * cell(ic)%ds

c...          Find data boundary values -- NEEDS WORK!:
              i2 = i0

              i3 = i1

c              WRITE(0,*) 'I2,3:',i2,i3,mode

c...          Flag...
c              tetarget = .FALSE.
c              IF ((index.EQ.1.OR.index.EQ.5).AND.osms28(i3,7).LT.0.0) 
c     .          tetarget = .TRUE.

              IF     (mode.EQ.1.OR.mode.EQ.2.OR.mode.EQ.3.OR.
     .                mode.EQ.5) THEN
c...            Interpolation boundary values provided in the input file:

c               *CRAP!*
                n0 = osmnode(i2)%ne  ! Needs work to allow relaxation 
                n1 = osmnode(i3)%ne
                p0 = osmnode(i2)%pe
                p1 = osmnode(i3)%pe
                t0 = osmnode(i2)%te
                t1 = osmnode(i3)%te

c                WRITE(0,*) 'N0,1:',n0,n1
c                WRITE(0,*) 'P0,1:',p0,p1
c                WRITE(0,*) 'T0,1:',t0,t1

                IF (osmnode(i2)%ne.EQ.-99.0.OR.
     .              osmnode(i2)%pe.EQ.-99.0.OR.
     .              osmnode(i2)%te.EQ.-99.0) THEN
c...              Linking to another plasma region where the solution has
c                 already (!) been calculated: 

                  CALL FindCell_New(i2,i3,it,icell)

c...              Assumptions: 1:1 mapping between cells and objects, the
c                 objects are 4 sided and have the standard DIVIMP indexing.  
c                 Also, cells in the core/SOL will always reference tubes
c                 that have a lower tube index via side 1-4 and PFZ cells 
c                 reference tubes with a higher index through side 2-3.
                  iobj = icell(2)
                  IF (tube(it)%type.EQ.GRD_PFZ) THEN
                    ic1 = obj(iobj)%omap(2)
                  ELSE
                    ic1 = obj(iobj)%omap(4)
                  ENDIF                  


c                  WRITE(0,*) ' MAP   ',iobj,ic1,nobj,ncell
                 
c...              Now have to search an see which tube this cell is in:
                  DO it1 = 1, ntube
                    IF (tube(it1)%cell_index(LO).LE.ic1.AND.
     .                  tube(it1)%cell_index(HI).GE.ic1) EXIT
                  ENDDO
                  IF (it1.EQ.ntube+1) 
     .              CALL ER('AssignNodeValues_New','Tube not '//
     .                      'identified',*99)



c                  WRITE(0,*) 'IC1,IT1',icell(2),cell(ic1)%ik
c                  WRITE(0,*) '       ',ic1,it1

                         
c                  IF (tube(it)%type.EQ.GRD_PFZ) THEN
c                    ik1 = ikouts(ikcell(2),ircell(2))  ! FIX
c                    ir1 = irouts(ikcell(2),ircell(2)) 
c                  ELSE
c                    ik1 = ikins(ikcell(2),ircell(2))   ! FIX
c                    ir1 = irins(ikcell(2),ircell(2))
c                  ENDIF


                  IF (osmnode(i2)%ne.EQ.-99.0) THEN
                    IF (density) THEN
                      n0 = fluid(ic1,1)%ne
                    ELSE
                      n0 = 2.0 * fluid(ic1,1)%te * fluid(ic1,1)%ne  ! Add Ti and M? 
                    ENDIF
                  ENDIF
                  IF (osmnode(i2)%pe.EQ.-99.0)   ! No v|| contribution!
     .              p0 = fluid(ic1,1)%ne * 
     .                   (fluid(ic1,1)%te + fluid(ic1,1)%ti)
                  IF (osmnode(i2)%te.EQ.-99.0) t0 = fluid(ic1,1)%te
c...              Shouldn't really be outer target (all this would go away if PSITARG was
c                 assigned properly):
                  IF (coord.EQ.3) psin0 = tube(it1)%psin
c                  IF (coord.EQ.4) rho0  = tube(it1)%rho

c                  WRITE(0,*) 't0:',t0
c                  WRITE(0,*) 'n0:',n0

                  WRITE(logfp,*) '  NODE TUBE LINK:',ic1,it1
                  WRITE(logfp,*) '            psin:',psin0
                  WRITE(logfp,*) '              ne:',n0  
                  WRITE(logfp,*) '              pe:',p0  
                  WRITE(logfp,*) '              te:',t0  
                ENDIF

c...            Make sure that t0,1 are positive:
                IF (tetarget) THEN
                  t0 = ABS(t0)
                  t1 = ABS(t1)
                ENDIF

                IF     (coord.EQ.1) THEN
c...              Linear along the line segment, nice and simple:
                  val0 = 0.0
                  val1 = 1.0
                  val = SNGL(tab)
                ELSEIF (coord.EQ.5) THEN
c...              Just give me PSIn:
                  psin0 = tube(it)%psin
                ELSE
c...            
                  IF ((osmnode(i2)%tube_range(1).NE.
     .                 osmnode(i3)%tube_range(1)).OR.
     .                (osmnode(i2)%tube_range(2).NE.
     .                 osmnode(i3)%tube_range(2)).OR.
     .                (osmnode(i2)%tube_range(1).GT.
     .                 osmnode(i2)%tube_range(2)))
     .              CALL ER('AssignNodeValues_2','Invalid '//
     .                      'tube index range, check input file',*99)

c...              Need range of PSIn over the segment:
                  DO it1 = osmnode(i2)%tube_range(1), 
     .                     osmnode(i2)%tube_range(2)
                    DO ic1 = tube(it1)%cell_index(LO), 
     .                       tube(it1)%cell_index(HI)
                      iobj = ic1
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
     .                    tcd.GE.0.0D0.AND.tcd.LT.1.0D0) THEN
                        IF     (tube(it)%type.EQ.GRD_SOL) THEN
                          IF (psin0.EQ.RHI) psin0 = tube(it1)%psin
                          IF (psin0.NE.RHI) psin1 = tube(it1)%psin
c                          WRITE(0,*) 'SOL:',psin0,psin1,it,it1
                        ELSEIF (tube(it)%type.EQ.GRD_PFZ) THEN
                          psin0 = MAX(psin0,tube(it1)%psin)
                          psin1 = MIN(psin1,tube(it1)%psin)
c                          WRITE(0,*) 'PFZ:',psin0,psin1,it,it1
                        ELSE
                          CALL ER('AssignNodeValues_2','Invalid '//
     .                            'RINGTYPE',*99)
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDDO

c                  WRITE(0,*) 'PSIN:',psin0,psin1,coord
                  IF     (coord.EQ.2) THEN
c...                Spatial along the IR-range of applicability:
                    STOP 'NOT READY 2'
                  ELSEIF (coord.EQ.3) THEN
c...                PSIn:
                    val0 = 0.0
                    val1 = ABS(psin1         - psin0)
                    val  = ABS(tube(it)%psin - psin0)
                    WRITE(logfp,*) ' psin:',psin0,psin1,tube(it)%psin
                  ELSEIF (coord.EQ.4) THEN
c...                RHO:
                    STOP 'NOT READY 4'
                  ELSE
                    CALL ER('AssignNodeValues_2','Invalid COORD A',*99)
                  ENDIF
                ENDIF

              ELSEIF (mode.EQ.4) THEN
c...            Load probe data, dummy values here:
                n0 = osmnode(i2)%ne
                p0 = osmnode(i2)%pe
                t0 = osmnode(i2)%te
                n1 = n0
                p1 = p0
                t1 = t0
              ELSE
                CALL ER('AssignNodeValues_2','Invalid MODE',*99)   
              ENDIF

c...          Check if quantities should be assigned:
              nc = .TRUE.
              pc = .TRUE.
              tc = .TRUE.
              IF (n0.EQ.0.0.OR.n1.EQ.0.0) nc = .FALSE.
              IF (p0.EQ.0.0.OR.p1.EQ.0.0) pc = .FALSE.
              IF (t0.EQ.0.0.OR.t1.EQ.0.0) tc = .FALSE.

              SELECTCASE (mode)
                CASE (1)
c...              Power law between v1 and v2:
                  val2 = (val - val0) / (val1 - val0)
                  IF (nc) ne(index) = n0 + val2**expon * (n1 - n0)
                  IF (pc) pe(index) = p0 + val2**expon * (p1 - p0)
                  IF (tc) te(index) = t0 + val2**expon * (t1 - t0)
                CASE (2)
c...              Exponential decay between v1 and v2:
                  C = expon  
                  IF (nc) THEN
                    A = (n1 - n0) / (EXP(-val1 / C) - 1.0)
                    B = n0 - A
                    ne(index) = A * EXP(-val / C) + B
                  ENDIF
                  IF (pc) THEN
                    A = (p1 - p0) / (EXP(-val1 / C) - 1.0)
                    B = p0 - A
                    pe(index) = A * EXP(-val / C) + B
                  ENDIF
                  IF (tc) THEN
                    A = (t1 - t0) / (EXP(-val1 / C) - 1.0)
                    B = t0 - A
                    te(index) = A * EXP(-val / C) + B
                  ENDIF
                CASE (3)
c...              Exponential decay to zero at infinity:
                  C = expon
                  A = n0 - n1
                  B = n1
                  IF (nc) ne(index) = A * EXP(-val / C) + B
                  A = p0 - p1
                  B = p1
                  IF (pc) pe(index) = A * EXP(-val / C) + B
                  A = t0 - t1
                  B = t1 
                  IF (tc) te(index) = A * EXP(-val / C) + B

                CASE (4)
c...              Load probe data from ascii file:
                  IF (nc) THEN
                    CALL LoadUpstreamData(osmnode(i1)%file_name,
     .                     itube,coord,osmnode(i1)%file_shift,
     .                     expon,osmnode(i2)%ne,tmp1) 
                    ne(index) = tmp1 * osmnode(i1)%file_scale_ne
                  ENDIF
                  IF (tc) THEN
                    CALL LoadUpstreamData(osmnode(i1)%file_name,
     .                     itube,coord,osmnode(i1)%file_shift,
     .                     expon,osmnode(i2)%te,tmp1) 
                    te(index) = tmp1 * osmnode(i1)%file_scale_te
                  ENDIF

                CASE (5)
c...              Interpolations from polynomial and exponential fitting parameters 
c                 that are listed in the input file:
 
c...              Count how many fit data lines are to be processed in this group:
                  nfit = 0
                  DO i4 = i1+1, osmnnode
                    IF (osmnode(i4)%type.NE.0.0) EXIT  
                    nfit = nfit + 1
                  ENDDO
c                  WRITE(logfp,*) 'NFIT:',nfit

                  DO i4 = 1, nfit
                    ifit = i1 + i4

c                    WRITE(logfp,*) 'FIT:',osmnode(ifit)%fit_type
c                    WRITE(logfp,*) '   :',osmnode(ifit)%fit_psin(1)
c                    WRITE(logfp,*) '   :',osmnode(ifit)%fit_psin(2)
c                    WRITE(logfp,*) '   :',osmnode(ifit)%fit_shift
c                    WRITE(logfp,*) '   :',osmnode(ifit)%fit_quantity
c                    WRITE(logfp,*) '   :',osmnode(ifit)%fit_p(1:6)

                    IF     (osmnode(ifit)%fit_type.EQ.-1.0) THEN 
                      psin1 = psin0 - osmnode(ifit)%fit_psin(1)
                    ELSEIF (psin1.GE.osmnode(ifit)%fit_psin(1).AND.
     .                      psin1.LE.osmnode(ifit)%fit_psin(2)) THEN
c...                  Only need this for the exponential fit coefficients
c                     returned in IDL (ts.pro at the moment):
                      IF (osmnode(ifit)%fit_type.EQ.2.0) THEN
                        psin2 = psin1 - osmnode(ifit)%fit_psin(1)
                      ELSE
                        psin2 = psin1
                      ENDIF
c                      WRITE(logfp,*) 'PSIN1,2:',psin1,psin2
                      SELECTCASE (NINT(osmnode(ifit)%fit_type))
                        CASE ( 1)  ! Polynomial
                          val = osmnode(ifit)%fit_p(1)            +
     .                          osmnode(ifit)%fit_p(2) * psin2    +
     .                          osmnode(ifit)%fit_p(3) * psin2**2 +
     .                          osmnode(ifit)%fit_p(4) * psin2**3 +
     .                          osmnode(ifit)%fit_p(5) * psin2**4 +
     .                          osmnode(ifit)%fit_p(6) * psin2**5
                        CASE ( 2)  ! Exponential
                          val = osmnode(ifit)%fit_p(1) * 
     .                          EXP(osmnode(ifit)%fit_p(2)*psin2) + 
     .                          osmnode(ifit)%fit_p(3)
                        CASE ( 3)  ! TANH 
                          p(0:4) = osmnode(ifit)%fit_p(1:5)
                          v0 = (p(0) - psin2) / (2.0 * p(1))               ! from /home/mastts/lib/edgefunctionats.pro
                          v1 = ((1. + p(3) * v0) * EXP(v0) - EXP(-v0)) /   ! from /home/mastts/bck/mtanh.pro (right one?) 
     .                                            (EXP(v0) + EXP(-v0))
                          val = (p(2) - p(4)) / 2.0 * (v1 + 1.0) + p(4)
                        CASEDEFAULT
                          CALL ER('AssignNodeValues_2',
     .                            'Unknown data type',*99)
                      ENDSELECT
                      SELECTCASE (NINT(osmnode(ifit)%fit_quantity))
                        CASE (1)  
                          IF (osmnode(ifit)%fit_type.EQ.3.0)  ! Special for TANH fit from ts.pro
     .                      val = val * 1.0E+19  
                          ne(index) = val
                        CASE (4)  
                          te(index) = val 
                        CASEDEFAULT
                          CALL ER('AssignNodeValues_2',
     .                            'Unknown data type',*99)
                      ENDSELECT
                    ENDIF
                  ENDDO
                CASE DEFAULT
                  CALL ER('AssignNodeValues_2','Invalid MODE',*99)   
              ENDSELECT

            ELSE
              CALL ER('AssignNodeValues_2','Invalid parameter '//
     .                'index',*99)
            ENDIF


c            IF (tetarget) THEN
c              te(index) = -te(index)
c            ENDIF


c...        Store node values:
            node_n = node_n + 1
            node_i(node_n) = index
            node_s(node_n)%s  = s (index)
            node_s(node_n)%ne = ne(index)
            node_s(node_n)%pe = pe(index)
            node_s(node_n)%te = te(index)

            node_s(node_n)%par_mode = osmnode(i3)%par_mode
            node_s(node_n)%par_exp  = osmnode(i3)%par_exp
            node_s(node_n)%par_set  = osmnode(i3)%par_set

          ENDIF

        ENDDO

      ENDDO



c... 

      node_n = node_n + 1
c...  Target nodes:
      node_s(1     )%s = 0.0
      node_s(node_n)%s = tube(it)%smax
      node_s(1     )%icell = 0
      node_s(node_n)%icell = tube(it)%n + 1
c...  Assign cell indeces:
      DO i1 = 2, node_n-1
        DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)
          IF (node_s(i1)%s.GE.cell(ic)%sbnd(1).AND.
     .        node_s(i1)%s.LE.cell(ic)%sbnd(2))THEN
            node_s(i1)%icell = ic - tube(it)%cell_index(LO) + 1
          ENDIF
        ENDDO
      ENDDO

      IF (log.GT.0) THEN
        WRITE(logfp,*) 
        DO i1 = 1, node_n
          WRITE(logfp,'(A,3I6,F10.2,2E10.2,F10.2)') 
     .      'NODE:',i1,node_i(i1),
     .      node_s(i1)%icell,node_s(i1)%s,
     .      node_s(i1)%ne,node_s(i1)%pe,node_s(i1)%te
        ENDDO
      ENDIF


c...  Check that only one symmetry node is specified:
      i2 = 0
      DO i1 = 1, node_n
        IF (node_i(i1).EQ.2) i2 = i2 + 1
      ENDDO

c...  Try to combine nodes:
      IF (i2.GT.1) THEN

c...    Find first instance:
        DO i1 = 1, node_n
          IF (node_i(i1).EQ.2) EXIT
        ENDDO          

c...    Check that all the nodes are at the same location:
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
c...        Check for quantity assigned multiple values:
            IF ((node_s(i1)%ne.NE.0.0.AND.node_s(i2)%ne.NE.0.0).OR.
     .          (node_s(i1)%pe.NE.0.0.AND.node_s(i2)%pe.NE.0.0).OR.
     .          (node_s(i1)%te.NE.0.0.AND.node_s(i2)%te.NE.0.0)) 
     .        CALL ER('AssignNodeValues_New','Degenerate '//
     .                'symmetry point node found',*99)
c...        Copy over data:
            IF (node_s(i2)%ne.NE.0.0) node_s(i1)%ne = node_s(i2)%ne
            IF (node_s(i2)%pe.NE.0.0) node_s(i1)%pe = node_s(i2)%pe
            IF (node_s(i2)%te.NE.0.0) node_s(i1)%te = node_s(i2)%te
c...        Delete degenerate node:
            IF (i2.LT.node_n) THEN
              node_s(i2:node_n-1) = node_s(i2+1:node_n)
              node_i(i2:node_n-1) = node_i(i2+1:node_n)
            ENDIF
            node_n = node_n - 1
          ENDIF
        ENDDO

      ENDIF



      ion = 1

c...  Check if target data is specified and overwrite existing tube 
c     data as required:
      DO i1 = 1, node_n 
        DO i2 = LO, HI
          IF (node_i(i1).EQ.2+i2) THEN
            IF (tarninter(i2).NE.0) THEN
              WRITE(0,*) 'WARNING AssignNodeValues: Overwriting '//
     .                   'specified target data'
              WRITE(0,*) '  TUBE = ',it,i2
            ENDIF
            IF (node_s(i1)%ne.LT.1.0E+10) THEN 
              tube(it)%jsat(i2,ion) = node_s(i1)%ne 
            ELSE
              tube(it)%ne(i2)     = node_s(i1)%ne
              tube(it)%ni(i2,ion) = node_s(i1)%ne 
            ENDIF
            tube(it)%te(i2)     = node_s(i1)%te
            tube(it)%ti(i2,ion) = node_s(i1)%te * opt%ti_ratio(i2)  ! *** Ti not loaded yet! ***
c...        Delete node from list:
            node_s(i1:node_n-1) = node_s(i1+1:node_n)
            node_i(i1:node_n-1) = node_i(i1+1:node_n)
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

      IF (log.GT.0) THEN
        DO i1 = 1, node_n
          WRITE(logfp,'(A,3I6,F10.2,2E10.2,F10.2)') 
     .      'NODE:',i1,node_i(i1),
     .      node_s(i1)%icell,node_s(i1)%s,
     .      node_s(i1)%ne,node_s(i1)%pe,node_s(i1)%te
        ENDDO
      ENDIF

c...  Set node indeces:
      DO i1 = 1, node_n 
        IF (node_i(i1).EQ.2) mnode = i1
      ENDDO
      nnode = node_n

c...  Assign values to nodes:
      node(1:nnode) = node_s(1:nnode)
c...  Assign other quantites:
      node(1:nnode)%jsat(1)   = 0.0
c      node(1:nnode)%pe        = 0.0
      node(1:nnode)%ni(1)     = 0.0
      node(1:nnode)%pi(1)     = 0.0
      node(1      :mnode)%ti(1) =node(1      :mnode)%te*opt%ti_ratio(LO) 
      node(mnode+1:nnode)%ti(1) =node(mnode+1:nnode)%te*opt%ti_ratio(HI) 
      node(1:nnode)%machno    = 0.0
      node(1:nnode)%potential = 0.0
      node(1:nnode)%efield    = 0.0

c      IF (node(1    )%ne   .EQ.0.) node(1    )%ne   =tube(it)%ne(LO)
c      IF (node(nnode)%ne   .EQ.0.) node(nnode)%ne   =tube(it)%ne(HI)
c      IF (node(1    )%ti(1).EQ.0.) node(1    )%ti(1)=tube(it)%ti(LO,1)
c      IF (node(nnode)%ti(1).EQ.0.) node(nnode)%ti(1)=tube(it)%ti(HI,1)


c...  Really don't like this but needed for C-Mod cases:
      DO i1 = 2, mnode
        IF (node(i1)%te.EQ.-98.0) node(i1)%te = node(i1-1)%te
      ENDDO
      DO i1 = mnode, nnode-1
        IF (node(i1)%te.EQ.-98.0) node(i1)%te = node(i1+1)%te
      ENDDO

c...  
      IF (tube(it)%te(LO).EQ.-99.0) THEN
        IF (node(2)%ne.EQ.0.0) 
     .    CALL ER('AssignNodeValues_New','Need density for sheath '//
     .            'limited particle flux calculation (LO)',*99)
        opt_tube%p_ion(LO) = 999
        tube(it)%te  (LO)     = node(2)%te
        tube(it)%ti  (LO,ion) = node(2)%ti(ion)
        tube(it)%jsat(LO,ion) = 
     .    GetJsat2(node(2)%te,node(2)%ti(ion),0.5*node(2)%ne,1.0) 
        STOP 'OK 1'
      ENDIF
c...  Set low index target node Te,i:
      SELECTCASE(node(2)%par_set)
        CASE (0) 
          node(1)%te = tube(it)%te(LO)
        CASE (1) 
          node(1)%te = node(2)%te
        CASE DEFAULT
          CALL ER('_New','Unknown LO PAR_SET',*99) 
      ENDSELECT
      node(1)%jsat(ion) = tube(it)%jsat(LO,ion)
      node(1)%ne        = 0.0
      node(1)%pe        = 0.0
      node(1)%ti(ion)   = tube(it)%ti(LO,ion)

      WRITE(logfp,*) 'JSAT9:',it,tube(it)%jsat(LO,ion)

      IF (tube(it)%te(HI).EQ.-99.0) THEN
        IF (node(nnode-1)%ne.EQ.0.0) 
     .    CALL ER('AssignNodeValues_New','Need density for sheath '//
     .            'limited particle flux calculation (HI)',*99)
        opt_tube%p_ion(HI) = 999
        tube(it)%te  (HI)     = node(nnode-1)%te
        tube(it)%ti  (HI,ion) = node(nnode-1)%ti(ion)
        tube(it)%jsat(HI,ion) = 
     .    GetJsat2(node(nnode-1)%te,node(nnode-1)%ti(ion),
     .             node(nnode-1)%ne*0.5,1.0) 
        STOP 'OK 2'
      ENDIF
      SELECTCASE(node(nnode-1)%par_set)
        CASE (0) 
          node(nnode)%te = tube(it)%te(HI)
        CASE (1) 
          node(nnode)%te = node(nnode-1)%te
        CASE DEFAULT
          CALL ER('_New','Unknown HI PAR_SET',*99) 
      ENDSELECT
      node(nnode)%jsat(ion) = tube(it)%jsat(HI,ion)
      node(nnode)%ne        = 0.0
      node(nnode)%pe        = 0.0
      node(nnode)%ti(ion)   = tube(it)%ti(HI,ion)
          
      IF (log.GT.0) THEN
        WRITE(logfp,*) 
        WRITE(logfp,*) 'NODE -:',nnode,mnode
        DO i1 = 1, node_n
          WRITE(logfp,'(A,4I6,F10.2,3E10.2,2F10.2)') 
     .      'NODE -:',i1,node_i(i1),node(i1)%par_set,
     .      node(i1)%icell,node_s(i1)%s,
     .      node(i1)%jsat(1),node(i1)%ne,
     .      node(i1)%pe,
     .      node(i1)%te,node(i1)%ti(1)
        ENDDO
      ENDIF

      RETURN
99    WRITE(0,*) ' ITUBE=',itube
      WRITE(0,*) ' MNODE=',mnode
      WRITE(0,*) ' NNODE=',nnode
      WRITE(0,*) ' TYPE =',tube(itube)%type
      WRITE(0,*) ' COORD=',coord
      STOP
      END


c
c ======================================================================
c
      SUBROUTINE AssignNodeValues_Legacy(itube,nnode,mnode,node)
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_geometry
      USE mod_legacy
      IMPLICIT none

      INTEGER itube,nnode,mnode       
      TYPE(type_node) :: node(*)

      REAL GetRelaxationFraction

      INTEGER ic,ic1,it,it1,i0,i1,i2,i3,i4,ifit,index,mode,
     .        iobj,isrf,ivtx(2)
      LOGICAL tc,nc,density,tetarget
      REAL    te(0:6),ne(0:6),s(0:6),
     .        frac,t0,t1,n0,n1,A,B,C,coord,expon,
     .        psin0,psin1,psin2,
     .        prb1,tmp1,val,val0,val1,val2,p(0:5),v0,v1,v2
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd


      INTEGER node_n,node_index(0:50)
      TYPE(type_node) :: node_store(0:50)


      s  = 0.0
      ne = 0.0
      te = 0.0

      node_n = 1

      frac = GetRelaxationFraction()

c...  Better/cleaner to pass the tube to this routine, and not need to
c     use mod_sol28_locals...?
      it = itube


      IF (osmns28.EQ.0) 
     .  CALL ER('AssignNodeValues_2','Profiles not found',*99)


      DO i1 = 2, osmns28
        i0 = i1 - 1
        IF ((osms28(i0,1).NE.osms28(i1,1)).OR.osms28(i1,1).EQ.0.0.OR.
     .      osms28(i1,2).EQ.0.0) CYCLE

c...    Do not apply data if IT is outside specified range:
        IF ((osms28(i1,11).NE.0.0.AND.REAL(it).LT.osms28(i1,11) ).OR.
     .      (osms28(i1,12).NE.0.0.AND.REAL(it).GT.osms28(i1,12)))
     .    CYCLE


c...    Need to keep track of iteration currently being assigned, so
c       that scans in upstream with step can be done serially, rather than
c       having everything on the same line...
c (Do you have any idea what I am talking about?)


c...    Check that rings from different grid regions are not in the same group
c       of rings:
        DO i2 = NINT(osms28(i1,11)), NINT(osms28(i1,12))-1
          IF (tube(i2)%type.NE.tube(i2+1)%type) THEN
            IF (log.GT.0.AND.tube(i2)%type.NE.GRD_CORE) THEN
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
c   5 - something strange...
c
c   coord = 1 - linear on line segment
c         = 2 - linear on line segment, but from first to last ring intersection
c         = 3 - PSIn 
c         = 4 - RHO
c

        index = NINT(osms28(i1,1)) 
        mode  = NINT(osms28(i1,2))
        coord = osms28(i1,3)
        expon = osms28(i1,4)

c...  Decide if specified upstream data is density or pressure:
        density = .TRUE.   ! *** I DON'T LIKE THIS SOLUTION ***
        IF (index.EQ.3.AND.
     .      (tube(it)%type.EQ.GRD_SOL.AND..FALSE..OR.
     .       tube(it)%type.EQ.GRD_PFZ.AND..FALSE.))
     .    density = .FALSE.


        a1 = DBLE(osms28(i1-1,5))
        a2 = DBLE(osms28(i1-1,6))
        b1 = DBLE(osms28(i1  ,5))
        b2 = DBLE(osms28(i1  ,6))

        DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)

c...      Assumed 1:1 mapping between grid and data:
          iobj = ic
          isrf = ABS(obj(iobj)%iside(1))
          ivtx(1:2) = srf(isrf)%ivtx(1:2)
          c1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
          c2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2)))
          isrf = ABS(obj(iobj)%iside(3))
          ivtx(1:2) = srf(isrf)%ivtx(1:2)
          d1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
          d2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2)))
c          id = korpg(ik,ir)
c          c1 = 0.5D0 * DBLE(rvertp(1,id) + rvertp(2,id))
c          c2 = 0.5D0 * DBLE(zvertp(1,id) + zvertp(2,id))
c          d1 = 0.5D0 * DBLE(rvertp(3,id) + rvertp(4,id))
c          d2 = 0.5D0 * DBLE(zvertp(3,id) + zvertp(4,id))
       
          CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
 
          IF (tab.GE.0.0D0.AND.tab.LT.1.0D0.AND.
     .        tcd.GE.0.0D0.AND.tcd.LT.1.0D0) THEN
c...        Intersecting between the line segment and the ring is found:

c...        These are here in case psin0 gets picked up when linking exponential
c           decay data to neighbouring ring:
            IF (tube(it)%type.EQ.GRD_PFZ) THEN
              psin0 = -RHI
              psin1 =  RHI
            ELSE
              psin0 =  RHI
              psin1 =  RHI
            ENDIF

            IF (.FALSE.) THEN
c...          Need a check to know if the ring is over-constrained:
            ELSEIF (index.GE.0.AND.index.LE.6) THEN

              IF (index.GE.1.AND.index.LE.5) THEN
                s(index) = cell(ic)%sbnd(1) + SNGL(tcd) * cell(ic)%ds
c                s(index)=ksb(ik-1,ir)+
c     .                   SNGL(tcd)*(ksb(ik,ir)-ksb(ik-1,ir))
              ENDIF

c...          Find data boundary values -- NEEDS WORK!:
              i2 = i0

              i3 = i1

c...          Flag...
              tetarget = .FALSE.
              IF ((index.EQ.1.OR.index.EQ.5).AND.osms28(i3,7).LT.0.0) 
     .          tetarget = .TRUE.

              IF     (mode.EQ.1.OR.mode.EQ.2.OR.mode.EQ.3.OR.
     .                mode.EQ.5) THEN
c...            Interpolation boundary values provided in the input file:

c               *CRAP!*
                t0 = osms28(i2,7) + frac*(osms28(i2,9) -osms28(i2,7))
                n0 = osms28(i2,8) + frac*(osms28(i2,10)-osms28(i2,8))
                t1 = osms28(i3,7) + frac*(osms28(i3,9) -osms28(i3,7))
                n1 = osms28(i3,8) + frac*(osms28(i3,10)-osms28(i3,8))

                IF ((index.EQ.3.OR.index.EQ.4.OR.index.EQ.5).AND.
     .              (osms28(i2,7).EQ.-99.0.OR.
     .               osms28(i2,8).EQ.-99.0)) THEN
c...              Linking to another plasma region where the solution has
c                 already (!) been calculated: 

                  STOP 'NEED NEW FINDCELL ROUTINE (PROPER IR REFERENCE)'
c                  CALL FindCell(i2,i3,ir,ikcell,ircell)
c                 ic1 = 
c                 it1 = 
                  IF (tube(it)%type.EQ.GRD_PFZ) THEN
c                    ik1 = ikouts(ikcell(2),ircell(2))  ! FIX
c                    ir1 = irouts(ikcell(2),ircell(2)) 
                  ELSE
c                    ik1 = ikins(ikcell(2),ircell(2))   ! FIX
c                    ir1 = irins(ikcell(2),ircell(2))
                  ENDIF


                  IF (osms28(i2,7).EQ.-99.0) t0 = fluid(ic1,1)%te
                  IF (osms28(i2,8).EQ.-99.0) THEN
                    IF (density) THEN
                      n0 = fluid(ic1,1)%ne
                    ELSE
                      n0 = 2.0 * fluid(ic1,1)%te * fluid(ic1,1)%ne  ! Add Ti and M? 
                    ENDIF
                  ENDIF
c...              Shouldn't really be outer target (all this would go away if PSITARG was
c                 assigned properly):
                  IF (coord.EQ.3) psin0 = tube(it1)%psin
c                  IF (coord.EQ.4) rho0  = tube(it1)%rho
                ENDIF

c...            Make sure that t0,1 are positive:
                IF (tetarget) THEN
                  t0 = ABS(t0)
                  t1 = ABS(t1)
                ENDIF

                IF     (coord.EQ.1) THEN
c...              Linear along the line segment, nice and simple:
                  val0 = 0.0
                  val1 = 1.0
                  val = SNGL(tab)
                ELSEIF (coord.EQ.5) THEN
c...              Just give me PSIn:
                  psin0 = tube(it)%psin
                ELSE
c...            
                  IF (NINT(osms28(i2,11)).NE.NINT(osms28(i3,11)).OR.
     .                NINT(osms28(i2,12)).NE.NINT(osms28(i3,12)).OR.
     .                NINT(osms28(i2,11)).GT.NINT(osms28(i2,12)))
     .              CALL ER('AssignNodeValues_2','Invalid '//
     .                      'IT range, check input file',*99)

c...              Need range of PSIn over the segment:
                  DO it1 = NINT(osms28(i2,11)), NINT(osms28(i2,12)) 
                    DO ic1 = tube(it1)%cell_index(LO), 
     .                       tube(it1)%cell_index(HI)
                      iobj = ic1
                      isrf = ABS(obj(iobj)%iside(1))
                      ivtx(1:2) = srf(isrf)%ivtx(1:2)
                      c1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
                      c2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2)))
                      isrf = ABS(obj(iobj)%iside(3))
                      ivtx(1:2) = srf(isrf)%ivtx(1:2)
                      d1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
                      d2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2)))
c                      id1 = korpg(ik1,ir1)
c                      c1 = 0.5D0 * DBLE(rvertp(1,id1) + rvertp(2,id1))
c                      c2 = 0.5D0 * DBLE(zvertp(1,id1) + zvertp(2,id1))
c                      d1 = 0.5D0 * DBLE(rvertp(3,id1) + rvertp(4,id1))
c                      d2 = 0.5D0 * DBLE(zvertp(3,id1) + zvertp(4,id1))
                      CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
                      IF (tab.GE.0.0D0.AND.tab.LT.1.0D0.AND.
     .                    tcd.GE.0.0D0.AND.tcd.LT.1.0D0) THEN
                        IF     (tube(it)%type.EQ.GRD_SOL) THEN
                          IF (psin0.EQ.HI) psin0 = tube(it1)%psin
                          IF (psin0.NE.HI) psin1 = tube(it1)%psin
                        ELSEIF (tube(it)%type.EQ.GRD_PFZ) THEN
                          psin0 = MAX(psin0,tube(it1)%psin)
                          psin1 = MIN(psin1,tube(it1)%psin)
                        ELSE
                          CALL ER('AssignNodeValues_2','Invalid '//
     .                            'RINGTYPE',*99)
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDDO

                  IF     (coord.EQ.2) THEN
c...                Spatial along the IR-range of applicability:
                    STOP 'NOT READY 2'
                  ELSEIF (coord.EQ.3) THEN
c...                PSIn:
                    val0 = 0.0
                    val1 = ABS(psin1         - psin0)
                    val  = ABS(tube(it)%psin - psin0)
                  ELSEIF (coord.EQ.4) THEN
c...                RHO:
                    STOP 'NOT READY 4'
                  ELSE
                    CALL ER('AssignNodeValues_2','Invalid COORD A',*99)
                  ENDIF
                ENDIF

              ELSEIF (mode.EQ.4) THEN
c...            Load probe data, dummy values here:
                t0 = osms28(i2,7)
                t1 = t0
                n0 = osms28(i2,8)
                n1 = n0
              ELSE
                CALL ER('S28params_v3','Invalid MODE',*99)   
              ENDIF

c...          Check if quantities should be assigned:
              tc = .TRUE.
              nc = .TRUE.
              IF (t0.EQ.0.0.OR.t1.EQ.0.0) tc = .FALSE.
              IF (n0.EQ.0.0.OR.n1.EQ.0.0) nc = .FALSE.

              SELECTCASE (mode)
                CASE (1)
c...              Power law between v1 and v2:
                  val2 = (val - val0) / (val1 - val0)
                  IF (tc) te(index) = t0 + val2**expon * (t1 - t0)
                  IF (nc) ne(index) = n0 + val2**expon * (n1 - n0)
                CASE (2)
c...              Exponential decay between v1 and v2:
                  C = expon  
                  IF (tc) THEN
                    A = (t1 - t0) / (EXP(-val1 / C) - 1.0)
                    B = t0 - A
                    te(index) = A * EXP(-val / C) + B
                  ENDIF
                  IF (nc) THEN
                    A = (n1 - n0) / (EXP(-val1 / C) - 1.0)
                    B = n0 - A
                    ne(index) = A * EXP(-val / C) + B
                  ENDIF
                CASE (3)
c...              Exponential decay to zero at infinity:
                  C = expon
                  A = t0 - t1
                  B = t1 
                  IF (tc) te(index) = A * EXP(-val / C) + B
                  A = n0 - n1
                  B = n1
                  IF (nc) ne(index) = A * EXP(-val / C) + B
                CASE (4)
c...              Load probe data from the .experiments file:
                  IF     (coord.EQ.3) THEN
                    prb1 = -1.0
                  ELSEIF (coord.EQ.4) THEN
                    prb1 = -3.0
                  ELSE
                    CALL ER('AssignNodeValues_2','Invalid COORD B',*99)   
                  ENDIF
                  IF (tc) THEN
                    STOP 'OPTION NOT READY A'
                    tmp1 = prb1
c                    CALL LoadProbeDataS28(ir,NINT(osms28(i2,7)),2,tmp1) 
                    te(index) = tmp1
                    IF     (expon.EQ.1.0) THEN
                      STOP 'OPTION NOT READY B'
                    ELSEIF (expon.EQ.2.0) THEN
                      te(0) = -te(index)
                      te(6) = -te(index)
                    ELSEIF (expon.EQ.3.0) THEN
                      te(0) = 97.0
                      te(6) = 97.0
                    ENDIF
                  ENDIF
                  IF (nc) THEN
                    STOP 'OPTION NOT READY C'
                    tmp1 = prb1
c                    CALL LoadProbeDataS28(ir,NINT(osms28(i2,8)),1,tmp1) 
                    ne(index) = tmp1
                    IF     (expon.EQ.1.0) THEN
                      STOP 'OPTION NOT READY D'
                    ELSEIF (expon.EQ.2.0) THEN
                      ne(0) = -ne(index)
                      ne(6) = -ne(index)
                    ELSEIF (expon.EQ.3.0) THEN
                    ENDIF
                  ENDIF

                CASE (5)
c...              Interpolations from polynomial and exponential fitting parameters 
c                 that are listed in the input file:
                  DO i4 = 1, NINT(osms28(i1,4))
                    ifit = i1 + i4
                    IF     (osms28(ifit,2).EQ.-1.0) THEN 
                      psin1 = psin0 - osms28(ifit,3)
                    ELSEIF (psin1.GE.osms28(ifit,3).AND.
     .                      psin1.LE.osms28(ifit,4)) THEN
c...                  Only need this for the exponential fit coefficients
c                     returned in IDL (ts.pro at the moment):
                      IF (osms28(ifit,2).EQ.2.0) THEN
                        psin2 = psin1 - osms28(ifit,3)
                      ELSE
                        psin2 = psin1
                      ENDIF
                       SELECTCASE (NINT(osms28(ifit,2)))
                        CASE ( 1)  ! Polynomial
                          val = osms28(ifit,7 )            +
     .                          osms28(ifit,8 ) * psin2    +
     .                          osms28(ifit,9 ) * psin2**2 +
     .                          osms28(ifit,10) * psin2**3 +
     .                          osms28(ifit,11) * psin2**4 +
     .                          osms28(ifit,12) * psin2**5
                        CASE ( 2)  ! Exponential
                          val = osms28(ifit,7) * 
     .                          EXP(osms28(ifit,8)*psin2) + 
     .                          osms28(ifit,9)
                        CASE ( 3)  ! TANH 
                          p(0:4) = osms28(ifit,7:11)
                          v0 = (p(0) - psin2) / (2.0 * p(1))               ! from /home/mastts/lib/edgefunctionats.pro
                          v1 = ((1. + p(3) * v0) * EXP(v0) - EXP(-v0)) /   ! from /home/mastts/bck/mtanh.pro (right one?) 
     .                                            (EXP(v0) + EXP(-v0))
                          val = (p(2) - p(4)) / 2.0 * (v1 + 1.0) + p(4)
                        CASEDEFAULT
                          CALL ER('AssignNodeValues_2',
     .                            'Unknown data type',*99)
                      ENDSELECT
                      SELECTCASE (NINT(osms28(ifit,6)))
                        CASE (1)  
                          IF (osms28(ifit,2).EQ.3.0) val = val * 1.0E+19  ! Special for TANH fit from ts.pro
                          ne(index) = val
                        CASE (4)  
                          te(index) = val 
                        CASEDEFAULT
                          CALL ER('AssignNodeValues_2',
     .                            'Unknown data type',*99)
                      ENDSELECT
                    ENDIF
                  ENDDO
                CASE DEFAULT
                  CALL ER('AssignNodeValues_2','Invalid MODE',*99)   
              ENDSELECT

            ELSE
              CALL ER('AssignNodeValues_2','Invalid parameter '//
     .                'index',*99)
            ENDIF


            IF (tetarget) THEN
              te(index) = -te(index)
            ENDIF

            IF (.NOT.density.AND.index.EQ.3.AND.nc.AND.mode.NE.4) THEN
c...          Convert the pressure value to density:
              IF (te(3).EQ.0.0) THEN
                CALL ER('AssignNodeValues_2','Te3 not assigned',*99)
              ELSE
c...            Assumes that the Mach no. is low (note: Ti.NE.Te is not a problem
c               since that situation is currently resolved in the SOL28 routine):
c                WRITE(0,*) 'CONVERTING PRESSURE TO DENSITY',ne(3)
c                WRITE(0,*) 'OLD D:',ne(3)
                ne(3) = ne(3) / (2.0 * te(3))
              ENDIF
            ENDIF

c...        Store node values:
            node_n = node_n + 1
            WRITE(0,*) 'INDEX:',index,s(index)
            node_index(node_n) = index
            node_store(node_n)%s  = s (index)
            node_store(node_n)%ne = ne(index)
            node_store(node_n)%te = te(index)

          ENDIF

        ENDDO

      ENDDO




      IF (.TRUE.) THEN
        node_n = node_n + 1
c...    Target nodes:
        node_store(1     )%s = 0.0
        node_store(node_n)%s = tube(it)%smax
c...    Assign cell indeces:
        DO i1 = 1, node_n
          DO ic = tube(it)%cell_index(LO),tube(it)%cell_index(HI)
            IF (node_store(i1)%s.GE.cell(ic)%sbnd(1).AND.
     .          node_store(i1)%s.LE.cell(ic)%sbnd(2))THEN
              node_store(i1)%icell = ic - tube(it)%cell_index(LO) + 1
c >>>> EH? >>>>              IF (i1.EQ.3) s(3) = cell(ic)%s
            ENDIF
          ENDDO
        ENDDO

        IF (log.GT.0) THEN
          WRITE(logfp,*) 
          DO i1 = 1, node_n
            WRITE(logfp,'(A,3I6,F10.2,E10.2,F10.2)') 
     .        'NODE:',i1,node_index(i1),
     .        node_store(i1)%icell,node_store(i1)%s,
     .        node_store(i1)%ne,node_store(i1)%te
          ENDDO
        ENDIF

c...    Check that only one symmetry node is specified:
        i2 = 0
        DO i1 = 1, node_n
          IF (node_index(i1).EQ.2) i2 = i2 + 1
        ENDDO
        IF (i2.NE.1) CALL ER('AssignNodeValues_Legacy','More than '//
     .                       'one symmetry point node found',*99)

c...    Sort nodes based on s-distance along the field line:
        DO i1 = 1, node_n-1
          DO i2 = i1+1, node_n
            IF (node_store(i2)%s.LT.node_store(i1)%s) THEN
              node_store(0)  = node_store(i1)
              node_index(0)  = node_index(i1)
              node_store(i1) = node_store(i2)
              node_index(i1) = node_index(i2)
              node_store(i2) = node_store(0)
              node_index(i2) = node_index(0)
            ENDIF
c...        Also check that there isn't more than one node in each cell:
            IF (node_store(i1)%icell.EQ.node_store(i2)%icell) 
     .        CALL ER('AssignNodeValues_Legacy','More than one node '//
     .                'detected in cell',*99)
          ENDDO
        ENDDO

        IF (log.GT.0) THEN
          WRITE(logfp,*) 
          DO i1 = 1, node_n
            WRITE(logfp,'(A,3I6,F10.2,E10.2,F10.2)') 
     .        'NODE:',i1,node_index(i1),
     .        node_store(i1)%icell,node_store(i1)%s,
     .        node_store(i1)%ne,node_store(i1)%te
          ENDDO
        ENDIF

c...    Set node indeces:
        DO i1 = 1, node_n 
          IF (node_index(i1).EQ.2) mnode = i1
        ENDDO
        nnode = node_n

c...    Assign values to nodes:
        node(1:nnode) = node_store(1:nnode)
c...    Assign other quantites:
        node(1:nnode)%jsat(1)   = 0.0
        node(1:nnode)%pe        = 0.0
        node(1:nnode)%ni(1)     = 0.0
        node(1:nnode)%pi(1)     = 0.0
        node(1:nnode)%ti(1)     = 0.0
        node(1:nnode)%machno    = 0.0
        node(1:nnode)%potential = 0.0
        node(1:nnode)%efield    = 0.0

        IF (node(1    )%ne   .EQ.0.) node(1    )%ne   =tube(it)%ne(LO)
        IF (node(nnode)%ne   .EQ.0.) node(nnode)%ne   =tube(it)%ne(HI)
        IF (node(1    )%te   .EQ.0.) node(1    )%te   =tube(it)%te(LO)
        IF (node(nnode)%te   .EQ.0.) node(nnode)%te   =tube(it)%te(HI)
        IF (node(1    )%ti(1).EQ.0.) node(1    )%ti(1)=tube(it)%ti(LO,1)
        IF (node(nnode)%ti(1).EQ.0.) node(nnode)%ti(1)=tube(it)%ti(HI,1)
c REAL FUNCTION GetJsat(te,ti,ne,v)
        node(1    )%jsat(1) = tube(it)%jsat(LO,1)
        node(nnode)%jsat(1) = tube(it)%jsat(HI,1)
        node(1    )%ne      = 0.0
        node(nnode)%ne      = 0.0

        WRITE(0,*) 'NODE:',nnode,mnode
        DO i1 = 1, node_n
          WRITE(0,'(A,3I6,F10.2,2E10.2,F10.2)') 
     .      'NODE:',i1,node_index(i1),
     .      node(i1)%icell,node_store(i1)%s,
     .      node(i1)%jsat(1),node(i1)%ne,node_store(i1)%te
        ENDDO

      ELSE

        s(0) = 0.0
        s(6) = tube(it)%smax
        DO i1 = 0, 6
          DO ic = tube(it)%cell_index(LO),tube(it)%cell_index(HI)
            IF (s(i1).GE.cell(ic)%sbnd(1).AND.
     .          s(i1).LE.cell(ic)%sbnd(2))THEN
              node(i1+1)%icell = ic - tube(it)%cell_index(LO) + 1
              IF (i1.EQ.3) s(3) = cell(ic)%s
            ENDIF
          ENDDO
        ENDDO

        nnode = 7
        mnode = 4

c...    Assign values to nodes:
        node(1:7)%s  = s (0:6)
        node(1:7)%ne = ne(0:6)
        node(1:7)%te = te(0:6)
c...    Assign other quantites:
        node(1:7)%jsat(1)   = 0.0
        node(1:7)%pe        = 0.0
        node(1:7)%ni(1)     = 0.0
        node(1:7)%pi(1)     = 0.0
        node(1:7)%ti(1)     = 0.0
        node(1:7)%machno    = 0.0
        node(1:7)%potential = 0.0
        node(1:7)%efield    = 0.0

        IF (node(1)%ne   .EQ.0.0) node(1)%ne    = tube(it)%ne(LO)
        IF (node(7)%ne   .EQ.0.0) node(7)%ne    = tube(it)%ne(HI)
        IF (node(1)%te   .EQ.0.0) node(1)%te    = tube(it)%te(LO)
        IF (node(7)%te   .EQ.0.0) node(7)%te    = tube(it)%te(HI)
        IF (node(1)%ti(1).EQ.0.0) node(1)%ti(1) = tube(it)%ti(LO,1)
        IF (node(7)%ti(1).EQ.0.0) node(7)%ti(1) = tube(it)%ti(HI,1)

c REAL FUNCTION GetJsat(te,ti,ne,v)
        node(1)%jsat(1) = tube(it)%jsat(LO,1)
        node(7)%jsat(1) = tube(it)%jsat(HI,1)
      ENDIF

      RETURN
99    CONTINUE
      STOP
      END



