c     -*-Fortran-*-
c
c ======================================================================
c
c OSM-EIRENE06 interface:
c
c Write input file
c Write triangle
c     -WriteEIRENE_06
c
c ======================================================================
c
      SUBROUTINE GetList(nlist,list,range)
      USE mod_sol28_global
      IMPLICIT none

      INTEGER  , INTENT(OUT) :: nlist,list(*)
      CHARACTER, INTENT(IN ) :: range*(*)

      INTEGER i,j,k,l,m,n,val1,val2,val3,o,r,s,istart,iend
      LOGICAL debug

      debug = .TRUE. 

      l = LEN_TRIM(range)

      nlist = 0

c...  Quick check to see if 'infinite range' has been set:
      IF (l.EQ.4.AND.range(1:4).EQ.'none') RETURN


      j = 0
      k = 0
      s = 0
      r = 0
      istart = 0
      iend   = 0
      DO i = 1, l
        IF (debug) WRITE(0,*) 'pass====>',i
        IF (range(i:i).EQ.'(') THEN
          j = i + 1
          istart = nlist + 1
          CYCLE
        ENDIF
        IF (range(i:i).EQ.')') k = i - 1
        IF (range(i:i).NE.','.AND.j.EQ.0) j = i
        IF (range(i:i).EQ.','.AND.j.NE.0) k = i - 1

        IF (k.GT.0.OR.i.EQ.l) THEN
          IF (i.EQ.1                      ) j = 1
          IF (i.EQ.l.AND.range(i:i).NE.')') k = l

          o = SCAN(range(j:k),"Rr")
          IF (o.NE.0) THEN
            IF (istart.EQ.0.OR.iend.EQ.0) 
     .        CALL ER('GetList','Brackets not set properly',*99)            
            o = o + j - 1
            IF (debug) WRITE(0,*) 'REPEAT:',range(o:k)
            READ(range(o+1:k),*) r     

            IF (debug) WRITE(0,*) 'REPEAT:',r,istart,iend,nlist
            
            DO o = 1, r
              list(nlist+1:nlist+1+(iend-istart+1)) = list(istart:iend)
              nlist = nlist + (iend-istart+1)
            ENDDO

            IF (debug) WRITE(0,*) 'REPEAT:',r,istart,iend,nlist
            istart = 0
            iend   = 0
          ELSE
            m = 0
            DO n = j+1, k
              IF (range(n:n).EQ.'-') m = n
            ENDDO
            IF (m.GT.0) THEN
              IF (debug) WRITE(0,*) 'J,M-1:',j,m-1
              IF (debug) WRITE(0,*) 'M+1,K:',m+1,k
              IF (debug) WRITE(0,*) 'VALS2: '//TRIM(range(j:m-1)) 
              IF (debug) WRITE(0,*) 'VALS2: '//TRIM(range(m+1:k))
              READ(range(j:m-1),*) val1
              READ(range(m+1:k),*) val2
              IF (val1.GT.val2) 
     .          CALL ER('GetList','Range not specified correctly',*99)            
              DO o = val1, val2
                nlist = nlist + 1
                list(nlist) = o
              ENDDO
              IF (debug) WRITE(0,*) 'VALS2 : ',nlist
            ELSE
              IF (debug) WRITE(0,*) 'J,K   :',j,k
              IF (debug) WRITE(0,*) 'VALS1 : '//TRIM(range(j:k))
              READ(range(j:k),*) val1
              nlist = nlist + 1
              list(nlist) = val1       
              IF (debug) WRITE(0,*) 'VALS1 : ',nlist
            ENDIF
          ENDIF
          j = 0
          k = 0
          s = 0
          r = 0
          IF (range(i:i).EQ.')') iend = nlist
        ENDIF
      ENDDO

      IF (debug) WRITE(0,*) 'list:',list(1:nlist)



      RETURN
 99   STOP
      END
c
c ======================================================================
c
      LOGICAL FUNCTION CheckIndex(index,range)
      USE mod_sol28_global
      IMPLICIT none

      INTEGER  , INTENT(IN) :: index
      CHARACTER, INTENT(IN) :: range*(*)

      INTEGER i,j,k,l,m,n,val1,val2,val3,o,r,s
      LOGICAL debug

      debug = .FALSE. 

      CheckIndex = .FALSE.

      l = LEN_TRIM(range)

c...  Quick check to see if 'infinite range' has been set:
      IF ((l.EQ.3.AND.range(1:3).EQ.'all')) THEN
c     .    (l.EQ.1.AND.range(1:1).EQ.'0'  ).OR.
c     .    (l.EQ.2.AND.range(1:2).EQ.'-1') THEN  - removed 22/02/2011, SL
        CheckIndex = .TRUE.
        RETURN
      ENDIF

c      IF (l.EQ.4.AND.range(1:4).EQ.'none') THEN
c        CheckIndex = .FALSE. 
c        RETURN
c      ENDIF


      IF (debug) WRITE(0,*) index,' >'//TRIM(range)//'< ',l
      j = 0
      k = 0
      s = 0
      r = 0
      DO i = 1, l
        IF (range(i:i).NE.','.AND.j.EQ.0) j = i
        IF (range(i:i).EQ.','.AND.j.NE.0) k = i - 1

        IF (k.GT.0.OR.i.EQ.l) THEN
          IF (i.EQ.1) j = 1
          IF (i.EQ.l) k = l

c          IF (debug) WRITE(0,*) TRIM(range(j:k))
c          IF (debug) WRITE(0,*) SCAN(range(j:k),"Rr")
          o = SCAN(range(j:k),"Rr")
          IF (o.NE.0) THEN
            o = o + j - 1
            IF (debug) WRITE(0,*) 'REPEAT:',range(o:k)
            READ(range(o+1:k),*) r     
            k = o - 1
            IF (debug) WRITE(0,*) 'REPEAT:',r
          ENDIF
          o = SCAN(range(j:k),"Ss")
          IF (o.NE.0) THEN
            o = o + j - 1
            IF (debug) WRITE(0,*) 'SKIP:',range(o:k)
            READ(range(o+1:k),*) s
            k = o - 1
            IF (debug) WRITE(0,*) 'SKIP:',s
          ENDIF

          m = 0
          DO n = j+1, k
c         DO n = j, k
            IF (range(n:n).EQ.'-') m = n
          ENDDO
          IF (m.GT.0) THEN
            IF (debug) WRITE(0,*) 'J,M-1:',j,m-1
            IF (debug) WRITE(0,*) 'M+1,K:',m+1,k
            IF (debug) WRITE(0,*) 'VALS2: '//TRIM(range(j:m-1)) 
            IF (debug) WRITE(0,*) 'VALS2: '//TRIM(range(m+1:k))
            READ(range(j:m-1),*) val1
            READ(range(m+1:k),*) val2
            DO o = 1, r+1
              IF (debug) WRITE(0,*) '     : ',val1,val2,index,o
              IF (index.GE.val1.AND.index.LE.val2) CheckIndex = .TRUE.
              val3 = val2 - val1
              val1 = val2 + 1 + s
              val2 = val1 + val3
            ENDDO
          ELSE
            IF (debug) WRITE(0,*) 'J,K   :',j,k
            IF (debug) WRITE(0,*) 'VALS1 : '//TRIM(range(j:k))
            READ(range(j:k),*) val1
            DO o = 1, r+1
              IF (debug) WRITE(0,*) '      : ',val1,index,o
              IF (index.EQ.val1) CheckIndex = .TRUE.
              val1 = val1 + 1 + s
            ENDDO
          ENDIF
          j = 0
          k = 0
          s = 0
          r = 0
        ENDIF
      ENDDO

      IF (debug) WRITE(0,*) 'CHECKINDEX:',CheckIndex
c      STOP 'test'
c      STOP


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE GetObjCentre(iobj,cen)
      USE mod_eirene06
c      USE mod_geometry
      IMPLICIT none

      INTEGER iobj,ivtx
c      INTEGER iobj,iside,isrf,ivtx
      REAL*8  cen(3),count


      cen = 0.0D0
      DO ivtx = 1, 3
        cen(1:3) = cen(1:3) + ver(tri(iobj)%ver(ivtx),1:3)
      ENDDO
      cen = cen / 3.0D0

c      count = 0.0D0
c      DO iside = 1, obj(iobj)%nside
c        isrf = ABS(obj(iobj)%iside(iside))
c        DO ivtx = 1, srf(isrf)%nvtx
c          count = count + 1.0D0
c          cen(1:3) = cen(1:3) + vtx(1:3,srf(isrf)%ivtx(ivtx))
c        ENDDO
c      ENDDO
c      cen = cen / count

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: WriteEireneObjects
c
      SUBROUTINE WriteEireneObjects
      USE mod_eirene06_parameters
      USE mod_eirene06
      USE mod_eirene06_locals
      USE mod_geometry
      USE mod_options
      IMPLICIT none

      REAL*8 gmCalcSurfaceArea  ! This will go in mod_geometry...

      INTEGER fp,ivtx,iobj
      LOGICAL found
      REAL    version

      INTEGER   ipts(4),i1,i2,v1,istart,iend,iside,isrf,iplasma,ibfield,
     .          iobj1(4),iside1(4),isrf1(4),
     .          ik,ir,it,ctardat,max_ik,max_ir,max_is,max_plasma,
     .          tar_maxik,tar_maxir,itarget
      LOGICAL   filter,warning
      REAL      tot_flux,frac,target_flux,frac1
      INTEGER, ALLOCATABLE :: tar_objects(:)
      REAL   , ALLOCATABLE :: tdata(:)      
      REAL*8 , ALLOCATABLE :: tar_area(:), tar_totarea(:,:)

      WRITE(eirfp,*) 'WRITING EIRENE OBJECT FILES'

      warning = .TRUE.

      version = 1.00

      fp = 99

      ALLOCATE(tdata(nobj))

      IF (photons.EQ.-1) THEN
c...    Load ionisation data from previous EIRENE call:
        CALL LoadTriangleData(7,0,13,0,tdata,'default')  
      ELSE
        tdata = -999.0
      ENDIF

c...  Dump vertices:
      OPEN(UNIT=fp,FILE='objects.npco_char',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      istart = 1
      iend   = nvtx
      IF (nvtxmap.NE.0) THEN
        WRITE(fp,*) nvtxmap
        DO ivtx = istart, iend
          IF (vtxmap(ivtx).NE.0) THEN
            WRITE(fp,'(I10,3F14.6)') vtxmap(ivtx),vtx(1:3,ivtx)*100.0D0 ! *** Increase accuracy...? ***
          ENDIF
        ENDDO
      ELSE 
        WRITE(fp,*) iend-istart+1
        DO ivtx = istart, iend
          WRITE(fp,'(I6,3F14.6)') ivtx-istart+1,vtx(1:3,ivtx)*100.0D0  ! *** Increase accuracy...? ***
        ENDDO
      ENDIF
      CLOSE(fp)      

c...  Dump sides:
      OPEN(UNIT=fp,FILE='objects.elemente',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      

      istart = 1
      iend   = nobj

      max_ik     = -1
      max_ir     = -1
      max_is     = -1
      max_plasma = -1
      DO iobj = istart, iend
        max_ik     = MAX(max_ik    ,obj(iobj)%index(IND_IK))
        max_ir     = MAX(max_ir    ,obj(iobj)%index(IND_IR))
        max_is     = MAX(max_is    ,obj(iobj)%index(IND_IS))
        max_plasma = MAX(max_plasma,obj(iobj)%index(IND_PLASMA))
      ENDDO

      WRITE(fp,'(5I9)') iend-istart+1,max_ik,max_ir,max_is,max_plasma
      DO iobj = istart, iend
c...    Get vertices associated with the base triangle (always side #1 if 
c       grid prepared properly):
        isrf = obj(iobj)%iside(1)
        IF (isrf.LT.0) THEN
c...      Need to reverse order of points:
          ipts(1) = srf(-isrf)%ivtx(3)  ! This must match the convention in CheckTetrahedronStructure
          ipts(2) = srf(-isrf)%ivtx(2)
          ipts(3) = srf(-isrf)%ivtx(1)
        ELSE
          ipts(1:3) = srf(isrf)%ivtx(1:3)
        ENDIF
c...    Forth point, or apex point of sorts, from any of the other sides:
        isrf = ABS(obj(iobj)%iside(2))
        DO ivtx = 1, 3
          IF (srf(isrf)%ivtx(ivtx).NE.ipts(1).AND.
     .        srf(isrf)%ivtx(ivtx).NE.ipts(2).AND.
     .        srf(isrf)%ivtx(ivtx).NE.ipts(3)) EXIT
        ENDDO
        IF (ivtx.EQ.4) 
     .    CALL ER('WriteEireneObjects','4th tetrahedron point '//
     .            'not identified',*99)
        ipts(4) = srf(isrf)%ivtx(ivtx)
        IF (nvtxmap.NE.0) THEN  
          DO ivtx = 1, 4
            IF (vtxmap(ipts(ivtx)).EQ.0) THEN
              STOP 'STOP: PROBLEM WITH VTXMAP'
            ELSE
              ipts(ivtx) = vtxmap(ipts(ivtx))
            ENDIF
          ENDDO          
        ENDIF
        WRITE(fp,'(I9,4X,4I8,4X,4I3,4X,4I6)') 
     .    iobj-istart+1,
     .    (ipts(v1),v1=1,4),  ! (tri(i1)%sideindex(1,v1),v1=1,3),
     .    (0       ,v1=1,4),  ! (tri(i1)%sideindex(2,v1),v1=1,3),
     .    obj(iobj)%index(IND_IK),
     .    obj(iobj)%index(IND_IR),
     .    obj(iobj)%index(IND_IS),
     .    obj(iobj)%index(IND_PLASMA)
      ENDDO
      CLOSE(fp)      

c...  Connection map:
      OPEN(UNIT=fp,FILE='objects.neighbors',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) nobj
      DO iobj = istart, iend
c...    Collect connection map information:
        DO iside = 1, obj(iobj)%nside
          iobj1 (iside) = obj(iobj)%omap(iside)                    
          iside1(iside) = obj(iobj)%smap(iside)
          isrf = ABS(obj(iobj)%iside(iside))
          isrf1 (iside) = srf(isrf)%index(IND_SURFACE)             ! Surface (block 2A)
        ENDDO
        IF (tetrahedrons.AND.  
     .      iobj1(1).EQ.0.AND.iside1(1).EQ.0.AND.isrf1(1).EQ.0) THEN   ! For toroidal boundary tetrahedral surfaces, which are "lost" at the moment...
          IF (surface(nsurface)%type   .EQ.NON_DEFAULT_STANDARD.AND.   ! Looking for the special tetrahedron catch all surface of desperation...
     .        surface(nsurface)%subtype.EQ.ADDITIONAL.AND.
     .        surface(nsurface)%index(1).EQ.-1) THEN
            isrf1(1) = surface(nsurface)%num
            IF (warning) THEN
              WRITE(0,*) 
              WRITE(0,*) '----------------------------------'
              WRITE(0,*) '   CLEANING TETRADRON MAP, BAD!'
              WRITE(0,*) '----------------------------------'
              WRITE(0,*) 
              warning = .FALSE.
            ENDIF
          ELSE
            STOP 'SORT OUT THIS ANNOYING BUSINESS, AGAIN...'
          ENDIF
c          isrf1(1) = 5  
c          STOP 'SORT OUT THIS ANNOYING BUSINESS'
        ENDIF
        IF (tetrahedrons.AND.  
     .      iobj1(2).EQ.0.AND.iside1(2).EQ.0.AND.isrf1(2).EQ.0) 
     .    WRITE(0,*) 'CRAP 2'
        IF (tetrahedrons.AND.  
     .      iobj1(3).EQ.0.AND.iside1(3).EQ.0.AND.isrf1(3).EQ.0) 
     .    WRITE(0,*) 'CRAP 3'
        IF (tetrahedrons.AND.  
     .      iobj1(4).EQ.0.AND.iside1(4).EQ.0.AND.isrf1(4).EQ.0) 
     .    WRITE(0,*) 'CRAP 4'

c        IF (tetrahedrons.AND.   ! *** THIS ONE WAS ACTIVE MOST RECENTLY *** -SL, 12.08.09
c    .       iobj1(1).EQ.0.AND.iside1(1).EQ.0.AND.isrf1(1).EQ.0) THEN   ! *HACK* temporary for MAST
c         isrf1(1) = 8  
c        ENDIF
c        IF (iobj1(2).EQ.0.AND.iside1(2).EQ.0.AND.isrf1(2).EQ.0) THEN   ! *HACK* temporary for toroidal surfaces...
c          isrf1(2) = 8 
c        ENDIF
c        IF (iobj1(3).EQ.0.AND.iside1(3).EQ.0.AND.isrf1(3).EQ.0) THEN   ! *HACK* temporary for toroidal surfaces...
c          isrf1(3) = 8 
c        ENDIF
c        IF (iobj1(4).EQ.0.AND.iside1(4).EQ.0.AND.isrf1(4).EQ.0) THEN   ! *HACK* temporary for toroidal surfaces...
c          isrf1(4) = 8 
c        ENDIF
c        WRITE(fp,'(I9,4X,4(I9,F14.3,I4,4X),2I6,2X,2I4)') iobj-istart+1,
        WRITE(fp,'(I9,4X,4(I9,I4,I4,4X),2I6,4X,2I4)') iobj-istart+1,
     .    (iobj1(v1),iside1(v1),isrf1(v1),v1=1,4),  
     .    obj(iobj)%index(IND_IK),
     .    obj(iobj)%index(IND_IR),
     .    0,0
c        WRITE(fp,'(13I10)') iobj-istart+1,
c     .    (iobj1(v1),iside1(v1),isrf1(v1),v1=1,4)
      ENDDO
      CLOSE(fp)      


      WRITE(eirfp,*) ' *** VFIELD AND BFIELD NEED CORRECTING ***'

c...  Plasma data:
      OPEN(UNIT=fp ,FILE='objects.plasma',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
c...  Header:
      WRITE(fp,'(A,F4.2,A)') 
     .  '* VERSION ',version,' OF '//
     .  fluid_code(1:LEN_TRIM(fluid_code))//
     .  ' PLASMA FILE FOR TETRAHEDRON GRID'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* BULK PLASMA DATA'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A9,2A8,2X,A10,2X,3A10,2X,3A12)') 
     .  '*   Index','Te','Ti','ne','vx',
     .  'vy','vz','Bx','By','Bz'
      WRITE(fp,'(A9,2A8,2X,A10,2X,3A10,2X,3A12)') 
     .  '*        ','(eV)','(eV)','(cm-3)','(cm s-1)',
     .  '(cm s-1)','(cm s-1)','(Tesla)','(Tesla)','(Tesla)'
      WRITE(fp,*) nobj
      DO iobj = istart, iend
        iplasma = obj(iobj)%index(IND_PLASMA)
        ibfield = obj(iobj)%index(IND_BFIELD)
        WRITE(fp,'(I9,2F8.2,2X,1P,E10.2,2X,3E10.2,2X,
     .             3E12.4,2X,E12.4,0P,6X,3I4,I6)') iobj-istart+1,
     .    plasma(1,iplasma),          ! Te (eV)
     .    plasma(2,iplasma),          ! Ti (eV)
     .    plasma(3,iplasma)*1.0E-06,  ! ne (cm-3)
     .    plasma(4,iplasma)*100.0,    ! vx (cm s-1)
     .    plasma(5,iplasma)*100.0,    ! vy 
     .    plasma(6,iplasma)*100.0,    ! vz
     .    bfield(1,ibfield),          ! Bx (Tesla)
     .    bfield(2,ibfield),          ! By
     .    bfield(3,ibfield),          ! Bz
     .    tdata(iobj),  ! Ionisation rate from previous run (w or w/o photons)...
     .    obj(iobj)%index(IND_IK),
     .    obj(iobj)%index(IND_IR),
     .    obj(iobj)%index(IND_IS),
     .    obj(iobj)%index(IND_PLASMA)
      ENDDO

      tot_flux = 0.0

c...  Target data:
      ctardat = 0
      tar_maxik = 0
      tar_maxir = 0
c...  First count up how many surfaces are targets:
      DO iobj = istart, iend
        IF (grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID) CYCLE
        DO iside = 1, obj(iobj)%nside  ! This will always be iside=1 for now...
          isrf = ABS(obj(iobj)%iside(iside))
          IF (srf(isrf)%index(IND_TARGET).NE.0) THEN 
            ctardat = ctardat + 1
            tar_maxik = MAX(tar_maxik,obj(iobj)%index(IND_IK))
            tar_maxir = MAX(tar_maxir,obj(iobj)%index(IND_IR))
          ENDIF
        ENDDO
      ENDDO
      IF (ctardat.EQ.0.OR.tar_maxik.EQ.0.OR.tar_maxir.EQ.0) THEN 
c...    So, no target segments were identified on the tetrahedral grid.  Check
c       if a region of interest was specified, and if so, assume that 
c       target recycling needs to be turned off:
        DO i1 = 1, opt_eir%tet_n
          IF (opt_eir%tet_type(i1).EQ.1.0) filter = .TRUE.
        ENDDO
        IF (.NOT.filter)
     .    CALL ER('WriteEireneObjects','Problem with target map',*99)
c       Seems that tetrahedron filters are active so assume that all target
c       segments have been excluded -- delete the target strata:
        WRITE(0,*) 
        WRITE(0,*) '--------------------------------------------------'
        WRITE(0,*) '         *** TURNING OFF TARGET STRATA ***        '
        WRITE(0,*) '--------------------------------------------------'
        WRITE(0,*) 
        DO i1 = nstrata, 1, -1
          IF (strata(i1)%distrib.EQ.'FFTFF') THEN
            DO i2 = i1, nstrata-1
              strata(i2) = strata(i2+1)
            ENDDO
            nstrata = nstrata - 1
          ENDIF
        ENDDO
      ENDIF        
c...  Dynamic allocation because this number could be large for tetrahedrons, and
c     then loop again over all objects, recording the relevant information:
      IF (ctardat.GT.0) THEN 
        ALLOCATE(tar_objects(ctardat))
        ALLOCATE(tar_area   (ctardat))
        ALLOCATE(tar_totarea(tar_maxik,tar_maxir))
        ctardat = 0
        tar_totarea = 0.0D0
        DO iobj = istart, iend
          IF (grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID) CYCLE
          DO iside = 1, obj(iobj)%nside 
            isrf = ABS(obj(iobj)%iside(iside))
            IF (srf(isrf)%index(IND_TARGET).NE.0) THEN
              ctardat = ctardat + 1
              tar_objects(ctardat) = iobj
              tar_area   (ctardat) = gmCalcSurfaceArea(isrf)
              ik = obj(iobj)%index(IND_IK)
              ir = obj(iobj)%index(IND_IR) 
              tar_totarea(ik,ir) =tar_totarea(ik,ir) + tar_area(ctardat)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      WRITE(fp,'(A)') '* TARGET DATA'
      WRITE(fp,*) ctardat
c helium
c      WRITE(fp,*) ctardat*2
c      DO iobj = istart, iend
c        IF (grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID) CYCLE
      DO itarget = 1, ctardat
        iobj = tar_objects(itarget)
        ik = obj(iobj)%index(IND_IK)
        ir = obj(iobj)%index(IND_IR) 
        DO iside = 1, obj(iobj)%nside
          isrf = ABS(obj(iobj)%iside(iside))
          IF (srf(isrf)%index(IND_TARGET).EQ.0) CYCLE
c...      Find corresponding target data, as set in ProcessFluidCode, based
c         on the fluid code cell/ring indices:
          found = .FALSE.
          DO it = 1, ntardat
            IF (tardat(it,2).NE.REAL(ik).OR.
     .          tardat(it,3).NE.REAL(ir)) CYCLE
            IF (.NOT.found) THEN
              found = .TRUE.         
              frac = SNGL(tar_area(itarget) / tar_totarea(ik,ir))

              iplasma = obj(iobj)%index(IND_PLASMA)
              ibfield = obj(iobj)%index(IND_BFIELD)  

c              WRITE(0,'(A,2I6,1P,7E10.2,0P,2F10.4)') 
c     .           'FLUX DATA:',ik,ir,
c     .            tardat(it,7 )*frac,tardat(it,13),
c     .            bfield(4,ibfield),
c     .            plasma(3,iplasma),tardat(it,5),tardat(it,10),
c     .           tardat(it,5)*tardat(it,10)*tardat(it,13)*
c     .           bfield(4,ibfield)*tar_area(itarget)*1.602E-19,
c     .           plasma(3,iplasma)/plasma(19,iplasma),
c     .           SQRT((plasma(1 ,iplasma)+plasma(2 ,iplasma))/
c     .                (plasma(17,iplasma)+plasma(18,iplasma)))

              IF (opt_fil%target_flux.EQ.1) THEN
                IF ( plasma(19,iplasma)                    .LT.0.001.OR.
     .              (plasma(17,iplasma)+plasma(18,iplasma)).LT.0.001)
     .            CALL ER('WriteEireneObjects','Error in target '//
     .                    'data scaling',*99)
                target_flux = tardat(it,5) * tardat(it,10) * 
     .                        tardat(it,13) *
     .                        bfield(4,ibfield) * tar_area(itarget) * 
     .                        1.6022E-19

                frac1 = plasma(3,iplasma)/plasma(19,iplasma) *
     .                  SQRT((plasma(1 ,iplasma)+plasma(2 ,iplasma))/
     .                       (plasma(17,iplasma)+plasma(18,iplasma)))

                target_flux = target_flux * (0.5 * (frac1 - 1.0) + 1.0)

                WRITE(88,'(A,2I6,1P,2E10.2,0P,2F10.4)') 
     .            'FLUX DATA:',ik,ir,
     .            tardat(it,7)*frac,target_flux,frac1,
     .            (0.5 * (frac1 - 1.0) + 1.0)
              ELSE
                target_flux = tardat(it,7)*frac
              ENDIF

c              WRITE(0,*) 'FRAC:',frac,ik,ir
              WRITE(fp,'(I9,I6,1P,E10.2,0P,2F8.2,1P,2E10.2,0P,
     .                   F6.2,1P,E10.2,0P,6X,3I4)') 
c...            Target quantities:
     .          iobj-istart+1,iside,    ! Triangle index and side index
     .          target_flux,            ! ion flux to surface for species 1 (Amps)
     .          tardat(it,6 ),          ! Te (eV)                                 
     .          tardat(it,8 ),          ! Ti (ev)
     .          tardat(it,9 )*1.0E-06,  ! ni (cm-3)
     .          tardat(it,10)*100.0,    ! v_para (cm s-1) (not read by EIRENE as yet)  
     .          tardat(it,11),          ! Mach no.        (not read)
     .          tardat(it,12)*1.0E-04,  ! jsat (A cm-2)   (not read)
     .          ik,ir,iside             ! Fluid grid indices, for debugging only

c helium
c              WRITE(fp,'(I9,I6,1P,E10.2,0P,2F8.2,1P,2E10.2,0P,
c     .                   F6.2,1P,E10.2,0P,6X,2I4)') 
c...            Target quantities:
c     .          iobj-istart+1,iside,    ! Triangle index and side index
c     .          tardat(it,7 ),          ! ion flux to surface for species 1 (Amps)
c     .          tardat(it,6 ),          ! Te (eV)
c     .          tardat(it,8 ),          ! Ti (ev)
c     .          tardat(it,9 )*1.0E-06,  ! ni (cm-3)
c     .          tardat(it,10)*100.0,    ! v_para (cm s-1) (not read by EIRENE as yet)  
c     .          tardat(it,11),          ! Mach no.        (not read)
c     .          tardat(it,12)*1.0E-04,  ! jsat (A cm-2)   (not read)
c     .          ik,ir                   ! Fluid grid indices, for debugging only
            ELSE
              CALL ER('WriteEireneObjects','Target data appears to  '//
     .                'be over-specified',*99)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      IF (ctardat.GT.0) THEN 
        DEALLOCATE(tar_objects)
        DEALLOCATE(tar_area)
        DEALLOCATE(tar_totarea)
      ENDIF
      CLOSE(fp)      

      WRITE(eirfp,*) 'DONE'

      IF (ALLOCATED(tdata )) DEALLOCATE(tdata )
      IF (ALLOCATED(vtxmap)) DEALLOCATE(vtxmap)
      IF (ALLOCATED(plasma)) DEALLOCATE(plasma)
      IF (ALLOCATED(bfield)) DEALLOCATE(bfield)

      RETURN
96    WRITE(0,*) 'WRITETRIANGEFILES: PROBLEMS WITH FILE ACCESS'
      STOP
99    WRITE(0,*) ' CTARDAT  =',ctardat
      WRITE(0,*) ' TAR_MAXIK=',tar_maxik
      WRITE(0,*) ' TAR_MAXIR=',tar_maxir
      STOP
      END
c
c
c ======================================================================
c
c  subroutine: BinItems
c
c
c ...don't do this on the fly because...
c
      SUBROUTINE BinItems(n,cen,nx,ny,nz,nlist,ilist,glist)
      USE mod_geometry
      IMPLICIT none

      INTEGER, INTENT(IN)  :: n,nx,ny,nz
      REAL                 :: cen(3,n)
      INTEGER, INTENT(OUT) :: nlist(nx,ny,nz),ilist(nx,ny,nz),glist(n)


      INTEGER i,ix,iy,iz,count
      REAL    xmin,xmax,ymin,ymax,zmin,zmax,x,z,dx,dy,dz

      REAL   , PARAMETER :: PI = 3.1415926536, TOL = 1.0E-07

c...  Convert x coordinate to r and z coordinate to phi:
      IF (.TRUE.) THEN
c     
c       z <---------|
c                 / |
c                /TH|
c               /   |
c              /   .|.
c             /     .
c                   x
c
         DO i = 1, n
           x = cen(1,i)
           z = cen(3,i)
           cen(1,i) = SQRT(x**2 + z**2)
           IF (ABS(x).LT.TOL) THEN
             IF (z.GT.0.0) cen(3,i) =       PI / 2.0
             IF (z.LT.0.0) cen(3,i) = 3.0 * PI / 2.0
           ELSE
             cen(3,i) = ATAN(z / x)
             IF (x       .LT.0.0) cen(3,i) =cen(3,i)+PI      ! Hopefully this is not
             IF (cen(3,i).LT.0.0) cen(3,i) =cen(3,i)+PI*2.0  ! compiler dependant...
           ENDIF
          cen(3,i) = cen(3,i) * 180.0 / PI
          IF (cen(3,i).GT.360.0+TOL) THEN
            CALL ER('BinItems','PHI > 360.0',*99)
          ENDIF
          IF ((ABS(cen(3,i)      ).LT.TOL).OR.
     .        (ABS(cen(3,i)-360.0).LT.TOL)) cen(3,i) = 0.0 
        ENDDO

      ENDIF


      nlist = 0
      ilist = 0
      glist = 0


      WRITE(geofp,*) '    BINNING'

c.... Find extents:
      xmin =  1.0E+20
      xmax = -1.0E+20
      ymin =  1.0E+20
      ymax = -1.0E+20
      zmin =  1.0E+20
      zmax = -1.0E+20
      DO i = 1, n
        xmin = MIN(xmin,cen(1,i))
        xmax = MAX(xmax,cen(1,i))
        ymin = MIN(ymin,cen(2,i))
        ymax = MAX(ymax,cen(2,i))
        zmin = MIN(zmin,cen(3,i))
        zmax = MAX(zmax,cen(3,i))
      ENDDO
      xmin = xmin * (1.0 - 0.001 * SIGN(1.0,xmin))
      xmax = xmax * (1.0 + 0.001 * SIGN(1.0,xmax))
      ymin = ymin * (1.0 - 0.001 * SIGN(1.0,ymin))
      ymax = ymax * (1.0 + 0.001 * SIGN(1.0,ymax))
      zmin = zmin * (1.0 - 0.001 * SIGN(1.0,zmin))
      zmax = zmax * (1.0 + 0.001 * SIGN(1.0,zmax))
      dx = (xmax - xmin) / REAL(nx)      
      dy = (ymax - ymin) / REAL(ny)      
      dz = (zmax - zmin) / REAL(nz)      
      IF (ABS(dz).LT.1.0E-6) dz = 1.0

c...  Assign zones (a bit painful doing it this way, but can't think
c     of how else to produce GLIST as a linear array, minimizing 
c     memory requirements):
      DO i = 1, n
        ix = INT((cen(1,i) - xmin) / dx) + 1
        iy = INT((cen(2,i) - ymin) / dy) + 1
        iz = INT((cen(3,i) - zmin) / dz) + 1
        nlist(ix,iy,iz) = nlist(ix,iy,iz) + 1         
      ENDDO

      count = 1
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            ilist(ix,iy,iz) = count
            count = count + nlist(ix,iy,iz)
          ENDDO
        ENDDO
      ENDDO

      WRITE(geofp,*) '    CHECK:',ilist(nx,ny,nz)+nlist(nx,ny,nz)-1,n

      nlist = 0
      DO i = 1, n
        ix = INT((cen(1,i) - xmin) / dx) + 1
        iy = INT((cen(2,i) - ymin) / dy) + 1
        iz = INT((cen(3,i) - zmin) / dz) + 1
        nlist(ix,iy,iz) = nlist(ix,iy,iz) + 1         
        glist(ilist(ix,iy,iz)+nlist(ix,iy,iz)-1) = i
      ENDDO      

      WRITE(geofp,*) '    DONE'


      RETURN
 99   STOP
      END
c
c ======================================================================
c
c  subroutine: RemoveDuplicateSurfaces
c
c
c ...don't do this on the fly because...
c
      SUBROUTINE RemoveDuplicateSurfaces
      USE mod_geometry
      IMPLICIT none

      INTEGER iobj,isrf,isrf1,ivtx,i1,i2,nx,ny,nz,ix,iy,iz,count,iside
      INTEGER, ALLOCATABLE :: nlist(:,:,:),ilist(:,:,:),glist(:),
     .                        mlist(:)

      REAL, ALLOCATABLE :: cen(:,:)  ! Does this really help memory management..?


      LOGICAL, PARAMETER :: aggressive = .FALSE.
      REAL   , PARAMETER :: PI = 3.1415926536, TOL = 1.0E-07


      WRITE(geofp,*) '  REMOVING DUPLICATE SURFACES'


! NEED SURFACE MAPPING INFO FIRST...  NO YOU DON'T SINCE NEEDED INFO
! IS RECORDED WHEN SURFACE IS CREATED -- OBJ ASSOCIATION! 

      SELECTCASE (0)
        CASE(0)

          WRITE(geofp,*) '    FINDING SURFACE CENTRES'

!...      Find geometric center of surface vertices:
          ALLOCATE(cen(3,nsrf))
          cen = 0.0D0
          DO isrf = 1, nsrf
            DO ivtx = 1, srf(isrf)%nvtx
              cen(1:3,isrf) =cen(1:3,isrf)+vtx(1:3,srf(isrf)%ivtx(ivtx))
            ENDDO
            cen(1:3,isrf) = cen(1:3,isrf) / REAL(srf(isrf)%nvtx)
          ENDDO


          nx = 30
          ny = 30
          nz = 23  ! 10

          ALLOCATE(nlist(nx,ny,nz))
          ALLOCATE(ilist(nx,ny,nz))
          ALLOCATE(glist(nsrf))
          ALLOCATE(mlist(nsrf))
          mlist = 0

          CALL BinItems(nsrf,cen,nx,ny,nz,nlist,ilist,glist)


          WRITE(geofp,*) '    SEARCHING FOR DUPLICATES'

c...      Search for matching surfaces:
          DO iz = 1, nz
            DO iy = 1, ny
              DO ix = 1, nx

                DO i1 = 0, nlist(ix,iy,iz)-2
                  isrf = glist(ilist(ix,iy,iz) + i1)

c                  IF (srf(isrf)%nvtx.NE.3) CYCLE
            
                  IF (mlist(isrf).NE.0) CYCLE
                  DO i2 = i1+1, nlist(ix,iy,iz)-1
                    isrf1 = glist(ilist(ix,iy,iz) + i2)
                    IF (MatchSurface(isrf,isrf1)) mlist(isrf1) = isrf
                  ENDDO
                ENDDO 

              ENDDO ! IX
            ENDDO   ! IY
          ENDDO     ! IZ

          ix = 0
          DO isrf = 1, nsrf
            IF (mlist(isrf).NE.0) ix = ix + 1
          ENDDO
          WRITE(geofp,*) '    ',ix,' DUPLICATES OF',nsrf,' FOUND'

          WRITE(geofp,*) '    UPDATING OBJECTS'

c...      Update object surface pointers:
          DO iobj = 1, nobj        
            DO iside = 1, obj(iobj)%nside            
              isrf = ABS(obj(iobj)%iside(iside))
              IF (mlist(isrf).NE.0) obj(iobj)%iside(iside)= -mlist(isrf)
            ENDDO
          ENDDO

c...      Update surface links:
          DO isrf = 1, nsrf
            IF (srf(isrf)%link.NE.0) THEN            
              WRITE(0,*) 'SOME WORK NEEDED HERE TO SORT OUT '//
     .                   'LINK UPDATES, STOPPING'
              STOP
            ENDIF
          ENDDO

c...      Delete extraneous surfaces:
          IF (aggressive) THEN
            WRITE(geofp,*) '    DELETING SURFACES'
          ENDIF

          IF (.FALSE.) THEN
            WRITE(geofp,*) '    CHECKING'
            DO isrf = 150000, nsrf-1
              DO isrf1 = isrf+1, nsrf
                IF (MatchSurface(isrf,isrf1).AND.mlist(isrf1).EQ.0) THEN
                  WRITE(geofp,*) 'CURSES: DUPLICATE FOUND'
                  WRITE(geofp,*) '  ISRF,1:',isrf,isrf1

                  DO iz = 1, nz
                    DO iy = 1, ny
                      DO ix = 1, nx
                        DO i1 = 0, nlist(ix,iy,iz)-1
                          IF (isrf .EQ.glist(ilist(ix,iy,iz)+i1)) 
     .                      WRITE(geofp,*) ' 0:',ix,iy,iz,cen(1:3,isrf)
                          IF (isrf1.EQ.glist(ilist(ix,iy,iz)+i1)) 
     .                      WRITE(geofp,*) ' 1:',ix,iy,iz,cen(1:3,isrf1)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO

                ENDIF
              ENDDO
            ENDDO
          ENDIF

        CASE DEFAULT
          STOP 'UNRECOGNIZED CASE IN DELETESURFACE'
      ENDSELECT

c...  Clear memory (move above once debugging is done?):
      IF (ALLOCATED(cen  )) DEALLOCATE(cen)     
      IF (ALLOCATED(nlist)) DEALLOCATE(nlist)
      IF (ALLOCATED(ilist)) DEALLOCATE(ilist)
      IF (ALLOCATED(glist)) DEALLOCATE(glist)
      IF (ALLOCATED(mlist)) DEALLOCATE(mlist)

      WRITE(geofp,*) '  DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c  subroutine: RemoveDuplicateVertices
c
      SUBROUTINE RemoveDuplicateVertices
      USE mod_geometry
      IMPLICIT none

      INTEGER ivtx,ivtx1,isrf,i1,count,nx,ny,nz,ix,iy,iz,i2,nremove
      INTEGER, ALLOCATABLE :: nlist(:,:,:),ilist(:,:,:),glist(:),
     .                        mlist(:)

      REAL, ALLOCATABLE :: cen(:,:)  ! Does this really help memory management..?

      LOGICAL, PARAMETER :: aggressive = .FALSE.

      WRITE(geofp,*) '  REMOVING DUPLICATE VERTICES'

!...  Find geometric center of surface vertices:
      ALLOCATE(cen(3,nvtx))
      cen = 0.0D0
      DO ivtx = 1, nvtx
        cen(1:3,ivtx) = vtx(1:3,ivtx)
      ENDDO


      nx = 30
      ny = 30
      nz = 23  ! 10

      ALLOCATE(nlist(nx,ny,nz))
      ALLOCATE(ilist(nx,ny,nz))
      ALLOCATE(glist(nvtx))
      ALLOCATE(mlist(nvtx))
      nlist = 0
      ilist = 0
      glist = 0
      mlist = 0

      CALL BinItems(nvtx,cen,nx,ny,nz,nlist,ilist,glist)

c...  Search for matching vertices:
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx

            DO i1 = 0, nlist(ix,iy,iz)-2
              ivtx = glist(ilist(ix,iy,iz) + i1)
              IF (mlist(ivtx).NE.0) CYCLE
              DO i2 = i1+1, nlist(ix,iy,iz)-1
                ivtx1 = glist(ilist(ix,iy,iz) + i2)

                IF (mlist(ivtx1).EQ.0.AND.
     .              DABS(vtx(1,ivtx)-vtx(1,ivtx1)).LT.1.0D-07.AND.
     .              DABS(vtx(2,ivtx)-vtx(2,ivtx1)).LT.1.0D-07.AND.
     .              DABS(vtx(3,ivtx)-vtx(3,ivtx1)).LT.1.0D-07) 
     .            mlist(ivtx1) = ivtx  ! Tag for reassignment

              ENDDO
            ENDDO 

          ENDDO ! IX
        ENDDO   ! IY
      ENDDO     ! IZ


      nremove = 0
      DO ivtx = 1, nvtx
        IF (mlist(ivtx).NE.0) nremove = nremove + 1
      ENDDO
      WRITE(geofp,*) '    REMOVING ',nremove,' OF',nvtx


c...  Replace reference to tagged vertices in surface lists:
      DO isrf = 1, nsrf
        DO ivtx = 1, srf(isrf)%nvtx            
          i1 = srf(isrf)%ivtx(ivtx)
          IF (mlist(i1).NE.0) srf(isrf)%ivtx(ivtx) = mlist(i1)
        ENDDO
        srf(isrf)%svtx = SUM(srf(isrf)%ivtx(1:srf(isrf)%nvtx))
      ENDDO

      nvtxmap = 0
      IF (aggressive) THEN
        WRITE(geofp,*) '    DELETING DUPLICATE VERTICES...'

c       Leaving this for now...
c        1-build list of vertecies to delete and pass to DeleteVertex
c        2-in DeleteVerted, build index map for non-deleted verteces
c        3-move data to get rid of delted vertices
c        4-do a single scan through surface definitions to get
c          rid of references to deleted verteces
c        the above bit of code can then be deleted I think...
c        adapt for deleteing surfaces quickly/efficiently... 
        DO ivtx = 1, nvtx
 
        ENDDO

c        count = 0
c        ivtx = nvtx + 1
c        DO WHILE (ivtx.GT.1)
c          ivtx = ivtx - 1
c          IF (mlist(ivtx).NE.0) THEN
c            count = count + 1
c            IF (MOD(count,1000).EQ.0) 
c     .        WRITE(geofp,*) '      ',count,' OF',nremove
c            CALL DeleteVertex(ivtx)
c            DO i1 = ivtx, nvtx-1
c              mlist(i1) = mlist(i1+1)
c            ENDDO
c          ENDIF
c        ENDDO
      ELSE
        IF (ALLOCATED(vtxmap)) DEALLOCATE(vtxmap)
        ALLOCATE(vtxmap(nvtx))
        vtxmap = 0
        DO ivtx = 1, nvtx
          IF (mlist(ivtx).EQ.0) THEN
            nvtxmap = nvtxmap + 1
            vtxmap(ivtx) = nvtxmap
          ENDIF 
        ENDDO
      ENDIF


c...  Clear memory (move above once debugging is done?):
      IF (ALLOCATED(cen)) DEALLOCATE(cen)     

      IF (ALLOCATED(nlist)) DEALLOCATE(nlist)
      IF (ALLOCATED(ilist)) DEALLOCATE(ilist)
      IF (ALLOCATED(glist)) DEALLOCATE(glist)
      IF (ALLOCATED(mlist)) DEALLOCATE(mlist)

      WRITE(geofp,*) '  DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c  subroutine: BuildConnectionMap
c
      SUBROUTINE BuildConnectionMap(istart,iend)
      USE mod_eirene06_locals
      USE mod_eirene06
      USE mod_geometry
      IMPLICIT none

      INTEGER, INTENT(IN) :: istart,iend

      INTEGER iobj,iside,isrf,ivtx,iobj1,iside1

      INTEGER, PARAMETER :: ns = 26
C     Krieger IPP/07 - SUN compiler chokes on this syntax, wants
C     variable declaration and initialization separately
      INTEGER s(3,ns)
C     INTEGER s(3,ns) / 0, 0,-1,  0,-1, 0, -1, 0, 0,  
      DATA    s       / 0, 0,-1,  0,-1, 0, -1, 0, 0,  
     .                  0, 0,+1,  0,+1, 0, +1, 0, 0, 
     .                 -1,-1, 0, +1,-1, 0, -1,+1, 0, +1,+1, 0, 
     .                 -1,-1,-1, +1,-1,-1, -1,+1,-1, +1,+1,-1, 
     .                 -1,-1,+1, +1,-1,+1, -1,+1,+1, +1,+1,+1, 
     .                 -1, 0,-1,  0,-1,-1, 
     .                 +1, 0,-1,  0,+1,-1,
     .                 -1, 0,+1,  0,-1,+1, 
     .                 +1, 0,+1,  0,+1,+1 /


      WRITE(geofp,*) '  BUILDING CONNECTION MAP'


c...  Clean up duplicate verticies, necessary for connection
c     map search:
      CALL RemoveDuplicateVertices
      CALL RemoveDuplicateSurfaces

c...  Removing duplicate surfaces essentially builds the
c     connection map since the entire grid is searched
c     to find matching surfaces, in order to eliminate
c     redundant data.  Note however that this will need to
c     be supplimented for local mesh refinement plasma
c     grids since neighbouring cells may no share common
c     surfaces (but all still fine for triangle/tetrahedron
c     grids):


      CALL BuildConnectionMap_New
      RETURN
 

c...  OLD CODE:
      WRITE(eirfp,*) '    MAPPING'
      DO iobj = istart, iend
        obj(iobj)%omap = 0
        obj(iobj)%smap = 0
      ENDDO
      DO iobj = istart, iend
        DO iside = 1, obj(iobj)%nside         
          isrf = obj(iobj)%iside(iside)
          IF (isrf.LT.0) THEN
            iobj1  = srf(-isrf)%obj
            iside1 = srf(-isrf)%side
            obj(iobj )%omap(iside ) = iobj1
            obj(iobj )%smap(iside ) = iside1
            obj(iobj1)%omap(iside1) = iobj 
            obj(iobj1)%smap(iside1) = iside
          ENDIF
        ENDDO
      ENDDO

      WRITE(eirfp,*) '  DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: BuildNewTriangleObjects
c
      SUBROUTINE BuildNewTriangleObjects
      USE mod_eirene06_parameters
      USE mod_eirene06
      USE mod_eirene06_locals
      USE mod_geometry
      IMPLICIT none

      TYPE(type_srf) newsrf

      INTEGER iver,itri,ivtx,isrf,iside
      REAL*8 a(3)

c...  Copy vertices:
      DO iver = 1, nver
        a(1:2) = ver(iver,1:2)
c        a(1:2) = DBLE(ver(iver,1:2))
        a(3)   = 0.0D0
        ntryvtx = ntryvtx + 1
        tryvtx(1:3,ntryvtx) = a(1:3)
c        ivtx = AddVertex(a)
      ENDDO 
c...  Build new triangle objects and surfaces:        

      
      ngrp = 2
      grp(1)%origin = GRP_MAGNETIC_GRID  ! *** NEED AddGroup and CopyGroup functions ***
      grp(1)%type   = GRP_TRIANGLE
      grp(2)%origin = GRP_VACUUM_GRID
      grp(2)%type   = GRP_TRIANGLE

      DO itri = 1, ntri
c...     
        try(itri)%group             = tri(itri)%type      ! Just works by luck...
        try(itri)%index             = 0
        try(itri)%index(IND_IK    ) = tri(itri)%index(1)
        try(itri)%index(IND_IR    ) = tri(itri)%index(2)
        try(itri)%index(IND_IS    ) = 0
        try(itri)%index(IND_ZONE  ) = tri(itri)%zone
c...    Plasma data:
        try(itri)%index(IND_PLASMA) = itri
        try(itri)%index(IND_BFIELD) = itri
        plasma(1 :20,itri) = tri(itri)%plasma(1:20)
        plasma(17:19,itri) = tri(itri)%plasma(1:3)   ! *** HACK *** target flux scaling for tetrahedrons...
        bfield(1:4 ,itri) = tri(itri)%bfield(1:4)
c...    Create surfaces:
        try(itri)%nside = 3

        DO iside = 1, try(itri)%nside
          try(itri)%omap(iside) = tri(itri)%map(iside)
          try(itri)%smap(iside) = tri(itri)%sid(iside)
c          try(itri)%map(iside) = REAL(tri(itri)%map(iside))+
c     .                           REAL(tri(itri)%sid(iside))/100.0
        ENDDO     

        newsrf%index(IND_STDGRD)  = tri(itri)%sideindex(1,1)
        newsrf%index(IND_TARGET)  = tri(itri)%sideindex(2,1)
        newsrf%index(IND_SURFACE) = tri(itri)%sur(1)
        newsrf%type    = SPR_LINE_SEGMENT
        newsrf%link    = 0
        newsrf%nvtx    = 2
        newsrf%ivtx(1) = tri(itri)%ver(1)
        newsrf%ivtx(2) = tri(itri)%ver(2)
        ntrysrf = ntrysrf + 1
        trysrf(ntrysrf) = newsrf
        try(itri)%iside(1) = ntrysrf

        newsrf%index(IND_STDGRD)  = tri(itri)%sideindex(1,2)
        newsrf%index(IND_TARGET)  = tri(itri)%sideindex(2,2)
        newsrf%index(IND_SURFACE) = tri(itri)%sur(2)
        newsrf%type    = SPR_LINE_SEGMENT
        newsrf%link    = 0
        newsrf%nvtx    = 2
        newsrf%ivtx(1) = tri(itri)%ver(2)
        newsrf%ivtx(2) = tri(itri)%ver(3)
        ntrysrf = ntrysrf + 1
        trysrf(ntrysrf) = newsrf
        try(itri)%iside(2) = ntrysrf

        newsrf%index(IND_STDGRD)  = tri(itri)%sideindex(1,3)
        newsrf%index(IND_TARGET)  = tri(itri)%sideindex(2,3)
        newsrf%index(IND_SURFACE) = tri(itri)%sur(3)
        newsrf%type    = SPR_LINE_SEGMENT
        newsrf%link    = 0
        newsrf%nvtx    = 2
        newsrf%ivtx(1) = tri(itri)%ver(3)
        newsrf%ivtx(2) = tri(itri)%ver(1)
        ntrysrf = ntrysrf + 1
        trysrf(ntrysrf) = newsrf
        try(itri)%iside(3) = ntrysrf
      ENDDO


      RETURN
 99   STOP
      END
c
c
c Read triangle
c Move data to regular grid
c
c ======================================================================
c
c subroutine: ProcessTetrahedrons_06
c
      SUBROUTINE ProcessTetrahedrons_06
      USE mod_eirene06_parameters
      USE mod_eirene06
      USE mod_eirene06_locals
      USE mod_geometry
      USE mod_options
      IMPLICIT none

      TYPE(type_srf   ) newsrf
      TYPE(type_object) newobj

      LOGICAL CheckIndex

c...  Bricks:

      INTEGER itry,nseg,i1,i2,isid,itet,ivtx,save_nobj,save_isrf,
     .        iobj,isrf,istart,iend,isector,isurface,idum1,
     .        isurface_list(nsurface)
      LOGICAL regular_grid
      LOGICAL :: debug  = .FALSE.
      LOGICAL :: hack   = .FALSE.
      LOGICAL :: filter = .FALSE.

      REAL*8 a(3,3),b(3,7),c(3,4),ang,ang1,ang2,dang(500,2),
     .       theta,frac,xcen,ycen,adelta

      REAL*8     DTOL
      PARAMETER (DTOL=1.0D-07)

C     Krieger IPP/07 - SUN compiler chokes on this syntax, wants
C     variable declaration and initialization separately
      INTEGER s(4,14)
C     INTEGER s(4,14) /3, 2, 1, 7,   
      DATA    s       /3, 2, 1, 7,   
     .                 1, 2, 5, 7,   5, 4, 1, 7,   
     .                 2, 3, 6, 7,   6, 5, 2, 7,   
     .                 3, 1, 4, 7,   4, 6, 3, 7,   
     .                 4, 5, 6, 7,
!...  Reversed:
     .                 1, 2, 4, 7,   4, 2, 5, 7,
     .                 2, 3, 5, 7,   5, 3, 6, 7,
     .                 3, 1, 6, 7,   6, 1, 4, 7 /

      INTEGER ishift,iside,itry1
      LOGICAL, ALLOCATABLE :: trycheck(:),trycycle(:)
      REAL   , ALLOCATABLE :: tryxcen(:) ,tryycen (:)

C     Krieger IPP/07 - SUN compiler chokes on this syntax, wants
C     variable declaration and initialization separately
      INTEGER t(3,4)
C     INTEGER t(3,4) /1, 2, 3,  1, 4, 2,  2, 4, 3,  3, 4, 1 /
      DATA    t      /1, 2, 3,  1, 4, 2,  2, 4, 3,  3, 4, 1 /

      WRITE(eirfp,*) 'BUILDING TETRAHEDRONS',debug,eirfp

      IF (.TRUE.) THEN
c...    Convert legacy triangle objects to generalized geometry objects:
        ntry = ntri
        ntrysrf = 0
        ntryvtx = 0
        ALLOCATE(try   (     ntry))
        ALLOCATE(trysrf(   3*ntry))
        ALLOCATE(tryvtx(3 ,6*ntry))
        ALLOCATE(plasma(20,  ntry))  ! Needs to be consistent with TMP_PLASMA in ResolveFilament
        ALLOCATE(bfield(4 ,  ntry))
        CALL BuildNewTriangleObjects
      ENDIF

      WRITE(eirfp,*) '  NTRY:',ntry
      WRITE(eirfp,*) '  NOBJ:',nobj
      WRITE(eirfp,*) '  NSRF:',nsrf
      WRITE(eirfp,*) '  NVTX:',nvtx

c...  Start of tetrahedron objects:
      istart = nobj + 1
      
c...  Build bricks:
      ngrp = 12
      grp(3)%origin = GRP_MAGNETIC_GRID  ! *** NEED AddGroup and CopyGroup functions ***
      grp(3)%type   = GRP_TETRAHEDRON
      grp(4)%origin = GRP_VACUUM_GRID
      grp(4)%type   = GRP_TETRAHEDRON

      regular_grid = .TRUE.

      IF (regular_grid) THEN   
c...    Toroidal distribution:
        ang = 360.0D0 / DBLE(ntorseg) * torfrac
        nseg = ntorseg
        dang(1,1) = 0.0D0
        dang(1,2) = dang(1,1) + ang
        DO isector = 2, nseg
          dang(isector,1) = dang(isector-1,2)
          dang(isector,2) = dang(isector  ,1) + ang
        ENDDO
      ELSE
        WRITE(0,*)
        WRITE(0,*) '=================================================='
        WRITE(0,*) '   WARNING: PARTICLE BALANCE ONLY SET FOR'
        WRITE(0,*) '            FOR REGULAR TOROIDAL GRID'
        WRITE(0,*) '=================================================='
        WRITE(0,*)
c        nseg = 12
c        dang(1 ,1) =  0.0D0    !  0.0D0       
c        dang(2 ,1) = 15.0D0    !  30.0D0      
c        dang(3 ,1) = 20.0D0    !  60.0D0      
c        dang(4 ,1) = 45.0D0    !  65.0D0      
c        dang(5 ,1) = 60.0D0    !  80.0D0      
c        dang(6 ,1) = 85.0D0    !  85.0D0      
c        dang(7 ,1) = 90.0D0    !  90.0D0      
c        dang(8 ,1) = 95.0D0    !  95.0D0      
c        dang(9 ,1) = 120.0D0   !  100.0D0     
c        dang(10,1) = 135.0D0   !  105.0D0     
c        dang(11,1) = 150.0D0   !  120.0D0     
c        dang(12,1) = 175.0D0   !  150.0D0     
c        dang(12,2) = 180.0D0   !  180.0D0     
        nseg = 12
        dang(1 ,1) = -30.0D0   !  0.0D0       
        dang(2 ,1) = -25.0D0   !  30.0D0      
        dang(3 ,1) = -20.0D0   !  60.0D0      
        dang(4 ,1) = -15.0D0   !  65.0D0      
        dang(5 ,1) = -10.0D0   !  80.0D0      
        dang(6 ,1) = -5.0D0    !  85.0D0      
        dang(7 ,1) =  0.0D0    !  90.0D0      
        dang(8 ,1) =  5.0D0    !  95.0D0      
        dang(9 ,1) =  10.0D0   !  100.0D0     
        dang(10,1) =  15.0D0   !  105.0D0     
        dang(11,1) =  20.0D0   !  120.0D0     
        dang(12,1) =  25.0D0   !  150.0D0     
        dang(12,2) =  30.0D0   !  180.0D0     
        DO isector = 1, nseg-1
          dang(isector,2) = dang(isector+1,1)
        ENDDO
      ENDIF

      adelta = dang(nseg,2) - dang(1,1)
      dang = dang - 0.5D0 * adelta      

      DO isector = 1, nseg
        WRITE(0,'(A,I6,2F12.6)') 
     .    'ANGLES:',isector,dang(isector,1:2)
        WRITE(eirfp,'(A,I6,2F12.6)') 
     .    'ANGLES:',isector,dang(isector,1:2)
      ENDDO

      dang = dang * D_DEGRAD

c      ang1 = 360.0D0 / DBLE(ntorseg) * D_DEGRAD

c...  Select which triangles will be part of the tetrahedron grid:
      ALLOCATE(trycycle(ntry))
      trycycle = .FALSE.
      DO itet = 1, opt_eir%tet_n
        IF (opt_eir%tet_type(itet).EQ.1.0) filter = .TRUE.
      ENDDO
      IF (filter) THEN
        trycycle = .TRUE.
        ALLOCATE(tryxcen(ntry))
        ALLOCATE(tryycen(ntry))
        DO itry = 1, ntry
          DO i1 = 1, 3
            isrf = try(itry)%iside(i1)
            i2 = 1
            IF (isrf.LT.0) i2 = 2 ! Side orientation is switched, so use other end point
            a(1,i1) = tryvtx(1,trysrf(ABS(isrf))%ivtx(i2))
            a(2,i1) = tryvtx(2,trysrf(ABS(isrf))%ivtx(i2))
            a(3,i1) = 0.0D0
          ENDDO
          tryxcen(itry) = SNGL(SUM(a(1,1:3))) / 3.0
          tryycen(itry) = SNGL(SUM(a(2,1:3))) / 3.0
        ENDDO
        DO itet = 1, opt_eir%tet_n
          IF (opt_eir%tet_type(itet).NE.1.0) CYCLE         
          WRITE(0,*) 'ITET:',itet
          WRITE(0,*) '    :',opt_eir%tet_type(itet)
          WRITE(0,*) '    :',opt_eir%tet_x1  (itet)
          WRITE(0,*) '    :',opt_eir%tet_y1  (itet)
          WRITE(0,*) '    :',opt_eir%tet_x2  (itet)
          WRITE(0,*) '    :',opt_eir%tet_y2  (itet)
          DO itry = 1, ntry
            IF (tryxcen(itry).GE.opt_eir%tet_x1(itet).AND.
     .          tryycen(itry).GE.opt_eir%tet_y1(itet).AND.
     .          tryxcen(itry).LE.opt_eir%tet_x2(itet).AND.
     .          tryycen(itry).LE.opt_eir%tet_y2(itet)) 
     .        trycycle(itry) = .FALSE.
          ENDDO
        ENDDO
        DEALLOCATE(tryxcen)
        DEALLOCATE(tryycen)
      ENDIF

      ALLOCATE(trycheck(ntry))
      trycheck = .FALSE.






      DO itry = 1, ntry  ! ntry   

        IF (trycycle(itry)) CYCLE

c        IF (try(itry)%index(IND_IK).NE.1.OR.
c     .      try(itry)%index(IND_IR).NE.5) CYCLE  ! separatrix on sonnet_13018_250.sm
c        IF (hack) CYCLE
c        hack = .TRUE.

c...    Assemble brick vertices:
        DO i1 = 1, 3
          isrf = try(itry)%iside(i1)
          i2 = 1
          IF (isrf.LT.0) i2 = 2 ! Side orientation is switched, so use other end point
          a(1,i1) = tryvtx(1,trysrf(ABS(isrf))%ivtx(i2))
          a(2,i1) = tryvtx(2,trysrf(ABS(isrf))%ivtx(i2))
          a(3,i1) = 0.0D0
        ENDDO

c...    Filter:
c        xcen = SUM(a(1,1:3)) / 3.0
c        ycen = SUM(a(2,1:3)) / 3.0
c        IF (xcen.LT.0.01) CYCLE
c        IF (xcen.LT.5.0.OR.ycen.LT.3.0) CYCLE
c        IF (ycen.GT.-2.0) CYCLE
c        IF (ycen.LT.3.0) CYCLE
c        IF (try(itry)%index(IND_IR).GT.1 .AND.
c     .      try(itry)%index(IND_IR).LT.20) CYCLE 
c        IF (ycen.LT.3.6.OR.xcen.LT.4.5) CYCLE        ! ITER 1
c        IF (ycen.LT.3.0) CYCLE
c        IF (ycen.LT.3.6.OR.xcen.LT.4.5) CYCLE        
c        IF (ycen.GT.0.2.OR.xcen.GT.0.7) CYCLE        ! C-Mod inner wall puff study


c...    Check orientation...?
        IF (debug) THEN
          WRITE(eirfp,*) 'A:',a(1:2,1)
          WRITE(eirfp,*) 'A:',a(1:2,2)
          WRITE(eirfp,*) 'A:',a(1:2,3)
        ENDIF

c...    Expand toroidally:
c        b(1,1:3) = a(1,1:3)
c        b(2,1:3) = a(2,1:3)
c        b(3,1:3) = a(1,1:3) * DTAN(-0.5D0*ang1)
c        b(1,4:6) = a(1,1:3)
c        b(2,4:6) = a(2,1:3)
c        b(3,4:6) = a(1,1:3) * DTAN(+0.5D0*ang1)

        isector = 1
        IF (debug) THEN
          WRITE(eirfp,*) 'DANG:',dang(isector,1:2)/D_DEGRAD
        ENDIF
        b(1,1:3) = a(1,1:3) * DCOS(dang(isector,1))
        b(2,1:3) = a(2,1:3)
        b(3,1:3) = a(1,1:3) * DSIN(dang(isector,1))
        b(1,4:6) = a(1,1:3) * DCOS(dang(isector,2))
        b(2,4:6) = a(2,1:3)
        b(3,4:6) = a(1,1:3) * DSIN(dang(isector,2))

c        b(1,1:3) = a(1,1:3) * DCOS(-dang(isector,1))
c        b(2,1:3) = a(2,1:3)
c        b(3,1:3) = a(1,1:3) * DSIN(-dang(isector,1))
c        b(1,4:6) = a(1,1:3) * DCOS(+dang(isector,2))
c        b(2,4:6) = a(2,1:3)
c        b(3,4:6) = a(1,1:3) * DSIN(+dang(isector,2))

c       Center of brick:
        b(1:3,7) = 0.0D0
        DO i1 = 1, 6
          b(1:3,7) = b(1:3,7) + b(1:3,i1) / 6.0D0
c        DO i1 = 1, 3                               ! Bug, sort of...
c          b(1:2,7) = b(1:2,7) + b(1:2,i1) / 3.0D0  
        ENDDO

c...    Check orientation...?
        IF (debug) THEN
          WRITE(eirfp,*) 'B:',b(1:3,1)
          WRITE(eirfp,*) 'B:',b(1:3,2)
          WRITE(eirfp,*) 'B:',b(1:3,3)
          WRITE(eirfp,*) 'B:',b(1:3,4)
          WRITE(eirfp,*) 'B:',b(1:3,5)
          WRITE(eirfp,*) 'B:',b(1:3,6)
          WRITE(eirfp,*) 'B:',b(1:3,7)
        ENDIF


c...    Determine group assignment (ad hoc at the moment):
c...    Make tetrahedrons:
        DO itet = 1, 8  ! 8 tetrahedrons for each triangle
          newobj = try(itry)
          newobj%segment(1) = itet
          newobj%group = try(itry)%group + 2 
          newobj%phi = SNGL(0.5D0*(dang(1,1) + dang(1,2)) / D_DEGRAD)
          newobj%nside = 4
          newobj%index(IND_IS) = 1
          newobj%index(IND_IK) = try(itry)%index(IND_IK)
          newobj%index(IND_IR) = try(itry)%index(IND_IR)

          ishift = 0
          IF (itet.GE.2.AND.itet.LE.7) THEN
            iside = INT(REAL(itet) / 2.0)
            itry1 = try(itry)%omap(iside)
c            itry1 = INT(try(itry)%map(iside))
            IF (itry1.NE.0) THEN
              IF (trycheck(itry1)) ishift = 7
            ENDIF
c            WRITE(eirfp,'(A,4I6,2I6)') 'TRYCHECK:',
c     .        itry,itet,iside,ishift,
c     .        try(itry)%smap(iside),itry1
          ENDIF

          c(1:3,1) = b(1:3,s(1,itet+ishift))
          c(1:3,2) = b(1:3,s(2,itet+ishift))
          c(1:3,3) = b(1:3,s(3,itet+ishift))
          c(1:3,4) = b(1:3,s(4,itet+ishift))

          IF (debug) THEN
            WRITE(eirfp,*)
            WRITE(eirfp,*) itet,ishift
            WRITE(eirfp,*) s(1:4,itet+ishift)
            WRITE(eirfp,*) 'C:',c(1:3,1)
            WRITE(eirfp,*) 'C:',c(1:3,2)
            WRITE(eirfp,*) 'C:',c(1:3,3)
            WRITE(eirfp,*) 'C:',c(1:3,4)
          ENDIF

c...      Assign index mapping:
          DO isid = 1, newobj%nside
c           Vertices:
            IF (isid.EQ.1.AND.itet.GE.2.AND.itet.LE.7) THEN
              i1 = INT(REAL(itet-1)/2.0+0.51)                
              IF (debug) WRITE(eirfp,*) 'IIII:',itet,i1
              newsrf%index = trysrf(ABS(try(itry)%iside(i1)))%index
          
              IF (debug) THEN
                WRITE(eirfp,*) ' :',
     .            trysrf(ABS(try(itry)%iside(i1)))%index(1:3)
              ENDIF
            ELSE
              newsrf%index = 0
            ENDIF
            newsrf%type = SPR_PLANAR_POLYGON
            newsrf%obj  = nobj + 1
            newsrf%side = isid
            newsrf%nvtx = 3
            DO ivtx = 1, 3
c              IF (DABS(c(1,t(ivtx,isid))).LT.1.0D-07) THEN
c                WRITE(0,*) 'TETRAHEDRONS: SUSPICIOUS X-VAL A',nobj+1
c                STOP 
c              ENDIF
              newsrf%ivtx(ivtx) = AddVertex(c(1,t(ivtx,isid)))
            ENDDO
            
            newobj%iside(isid) = AddSurface(newsrf)

            IF (debug) THEN
              WRITE(eirfp,*) 'I.:',isid
              WRITE(eirfp,*) 'IN:',newsrf%index(1:3)
              WRITE(eirfp,*) 'I0:',t(1:3,isid)
              WRITE(eirfp,*) 'I1:',newsrf%ivtx(1:3)
              WRITE(eirfp,*) 'I1:',srf(ABS(newobj%iside(isid)))%
     .                                     ivtx(1:3)

              WRITE(eirfp,*) 'V:',vtx(1:2,newsrf%ivtx(1))
              WRITE(eirfp,*) 'V:',vtx(1:2,newsrf%ivtx(2))
              WRITE(eirfp,*) 'V:',vtx(1:2,newsrf%ivtx(3))
            ENDIF

          ENDDO
c...      Center of tetrahedron -- needs to be done properly:
          newobj%x = SNGL(0.25D0 * SUM(c(1,1:4)))
          newobj%y = SNGL(0.25D0 * SUM(c(2,1:4)))
          newobj%z = SNGL(0.25D0 * SUM(c(3,1:4)))

c          IF (itry.EQ.1) THEN
c            WRITE(0,*) 'Cx:',c(1,1:4)
c            WRITE(0,*) 'Cy:',c(2,1:4)
c            WRITE(0,*) 'Cz:',c(3,1:4)
c            WRITE(0,*) 'Cen:',newobj%x,newobj%y,newobj%z
c          ENDIF

c...      Add object:
          idum1 = AddObject(newobj)

        ENDDO  ! Tetrahedron loop

        trycheck(itry) = .TRUE.
         
      ENDDO  ! Triangle loop

c...  Done with these, clear some memory:
      DEALLOCATE(trycheck)
      DEALLOCATE(trycycle)
      DEALLOCATE(try)
      DEALLOCATE(trysrf)
      DEALLOCATE(tryvtx)

      WRITE(eirfp,*) '  NTRY:',ntry
      WRITE(eirfp,*) '  NOBJ:',nobj
      WRITE(eirfp,*) '  NSRF:',nsrf
      WRITE(eirfp,*) '  NVTX:',nvtx
      WRITE(eirfp,*) '  TOROIDAL REPLICATION'

c...  Toroidal replication:
      save_nobj = nobj
c      DO ang = dang,  0.0D0*D_DEGRAD, dang
c      DO ang = 179.9D0*D_DEGRAD, 359.9D0*D_DEGRAD, dang
c      DO ang = dang, 179.9D0*D_DEGRAD, dang
c      DO ang = dang, 359.9D0*D_DEGRAD, dang
c      DO ang = DBLE(torus1)+dang, DBLE(torus2*0.99999)*D_DEGRAD, ang1

c      isector = 1   
c      DO ang = DBLE(torus1)+ang1, DBLE(torus2*0.99999)*D_DEGRAD, ang1
c        isector = isector + 1

      DO isector = 2, nseg ! nseg  
        ang1 = dang(isector,1) - dang(1,1)
        ang2 = dang(isector,2) - dang(1,2)
c        ang  = 0.5D0 * (ang1+ ang2) 
c        ang1 = ang - (dang(isector,1) - dang(1,1)) 
c        ang2 = ang + (dang(isector,2) - dang(1,2)) 

        WRITE(eirfp,'(A,I6,3F12.6)') '   TOROIDALIZING:',isector,
     .    SNGL(0.5D0 * (dang(isector,1)+dang(isector,2)) / D_DEGRAD),
     .    ang1/D_DEGRAD,ang2/D_DEGRAD

        DO iobj = 1, save_nobj ! save_nobj
          newobj = obj(iobj)
          newobj%index(IND_IS) = isector
          newobj%phi = 
     .      SNGL(0.5D0*(dang(isector,1) + dang(isector,2)) / D_DEGRAD)
c          newobj%phi = SNGL(ang) / D_DEGRAD 

          DO isid = 1, newobj%nside
            isrf = ABS(newobj%iside(isid))  ! I don't need to do any trickery here to make sure
            save_isrf = 0                   ! that new sides are CCW from the outside of the object
                                            ! because the properly oriented side will be duplicated
            DO WHILE (.TRUE.)               ! first and the subsequent CW sides will then be mapped 
                                            ! to those sides...
              newsrf = srf(isrf)
              newsrf%obj  = nobj + 1
              newsrf%side = isid
              DO i1 = 1, newsrf%nvtx
c...            Find toroidal angle of vertex relative to the minimum
c               toroidal angle in sector 1, since these tetrahedrons
c               may be stretched

                a(1:3,1) = vtx(1:3,newsrf%ivtx(i1))

                IF (DABS(a(1,1)).LT.DTOL) THEN
c                IF (a(1,1).EQ.0.0D0) THEN
                  ang = 0.0D0
                ELSE
                  ang = DATAN(a(3,1) / a(1,1))
                ENDIF
          
c                WRITE(0,'(A,I9,I4,12X,4F12.4)') 
c     .            'ANG:',iobj,i1,ang/D_DEGRAD,a(1:3,1)

                frac = (ang - dang(1,1)) / (dang(1,2) - dang(1,1))

                ang = (1.0D0 - frac) * (dang(isector,1) - dang(1,1)) + 
     .                         frac  * (dang(isector,2) - dang(1,2))

c                WRITE(0,'(A,I9,I4,2F12.4)') 
c     .            '   :',iobj,i1,frac,ang/D_DEGRAD

                a(1,2) = DCOS(ang) * a(1,1) - DSIN(ang) * a(3,1)
                a(2,2) = a(2,1)
                a(3,2) = DSIN(ang) * a(1,1) + DCOS(ang) * a(3,1)
c                IF (DABS(a(1,2)).LT.1.0D-07) THEN
c                  WRITE(0,*) 'TETRAHEDRONS: SUSPICIOUS X-VAL B',nobj+1
c                  STOP 
c                ENDIF
                newsrf%ivtx(i1) = AddVertex(a(1,2))

c                WRITE(0,'(A,I9,I4,24X,3F12.4)') 
c     .            '   :',iobj,i1,a(1:3,1) 

c                a(1:3,1) = vtx(1:3,newsrf%ivtx(i1))
c                a(1  ,2) = DCOS(ang) * a(1,1) - DSIN(ang) * a(3,1)
c                a(2  ,2) = a(2,1)
c                a(3  ,2) = DSIN(ang) * a(1,1) + DCOS(ang) * a(3,1)
c                newsrf%ivtx(i1) = AddVertex(a(1,2))
              ENDDO


              IF (save_isrf.EQ.0) THEN 
                newobj%iside(isid) = AddSurface(newsrf)
              ELSE
                isrf = AddSurface(newsrf)
                srf(save_isrf)%link = isrf
              ENDIF

              IF (srf(isrf)%link.EQ.0) THEN
                EXIT
              ELSE
                save_isrf = isrf
                isrf = srf(isrf)%link 
                WRITE(eirfp,*) 'ACTUALLY USING THIS CODE?'
                STOP
              ENDIF

            ENDDO

          ENDDO

c...      Barf, need to tidy this up since it's a repeat of the above code:
          a(1,1) = DBLE(obj(iobj)%x)
          a(2,1) = DBLE(obj(iobj)%y)
          a(3,1) = DBLE(obj(iobj)%z)
          IF (DABS(a(1,1)).LT.DTOL) THEN
            ang = 0.0D0
          ELSE
            ang = DATAN(a(3,1)/a(1,1))
          ENDIF
          frac = (ang - dang(1,1)) / (dang(1,2) - dang(1,1))
          ang = (1.0D0 - frac) * (dang(isector,1) - dang(1,1)) + 
     .                   frac  * (dang(isector,2) - dang(1,2))
          a(1,2) = DCOS(ang) * a(1,1) - DSIN(ang) * a(3,1)
          a(2,2) = a(2,1)
          a(3,2) = DSIN(ang) * a(1,1) + DCOS(ang) * a(3,1)
          newobj%x = SNGL(a(1,2))
          newobj%y = SNGL(a(2,2))
          newobj%z = SNGL(a(3,2))
          IF (iobj.EQ.1) THEN
            WRITE(0,*) 'Cen:',newobj%x,newobj%y,newobj%z
          ENDIF

          idum1 = AddObject(newobj)

c          if (iobj.Eq.10)    STOP 'sdfsdfd'

        ENDDO
      ENDDO

 

C     IN RAY/OUT, I DON'T YET HAVE A SITUATION WHERE ORIENTATION IS A BIG ISSUE, SO FOR EIRENE
c     I'LL JUST DO THE ORDERING ON THE FLY, CHECKING IF THE POINTS ARE CO- OR COUNTER CLOCKWISE
c     WITH RESPECT TO THE CENTER OF THE TETRAHEDRON, BUT NEED TO COME UP WITH A TEST FIRST, 
c     PERHAPS THE NORMAL WITH RESPECT TO THE CENTER VECTOR? 

      iend = nobj

      WRITE(eirfp,*) '  NOBJ:',nobj
      WRITE(eirfp,*) '  NSRF:',nsrf
      WRITE(eirfp,*) '  NVTX:',nvtx

c...  Build connection map:
      CALL BuildConnectionMap(istart,iend)

c.... For Eirene, need to make sure the tetrahedrons are arranged according to
c     convention:
      CALL FixTetrahedrons(istart,iend)

c...  Plasma association:
c      CALL CheckTetrahedronStructure  ! Add this here..?

c...  Add slices here, if requested:



c...  Impose filament structures:
      IF (opt_fil%opt.NE.0) THEN
        WRITE(eirfp,*) '  NOBJ:',nobj
        WRITE(eirfp,*) '  NSRF:',nsrf
        WRITE(eirfp,*) '  NVTX:',nvtx
        CALL ResolveFilament(-1)
        CALL AssignFilamentPlasma   
        iend = nobj
        CALL BuildConnectionMap(istart,iend)
        CALL FixTetrahedrons   (istart,iend)
      ENDIF

c...  Manual refinement:
      DO itet = 1, opt_eir%tet_n
        IF (opt_eir%tet_type(itet).NE.5.0) CYCLE
        CALL ResolveFilament(itet)
        iend = nobj
        CALL BuildConnectionMap(istart,iend)
        CALL FixTetrahedrons   (istart,iend)
      ENDDO

      CALL CheckTetrahedronStructure

      IF (.TRUE.) THEN
c...    Make a linked list for the non-default standard surfaces in the surface 
c       array (which also includes wall surfaces, which are not of interest here):
        isurface_list = 0
        DO i1 = 1, nsurface
          IF (surface(i1)%type.NE.NON_DEFAULT_STANDARD) CYCLE
          isurface_list(surface(i1)%num) = i1
        ENDDO

c        WRITE(0,*) 'LIST:'
c        WRITE(0,*) isurface_list(1:nsurface)

        DO iobj = 1, nobj
c...      Collect connection map information:
          DO iside = 1, obj(iobj)%nside
            isrf   = ABS(obj(iobj)%iside(iside))

            isurface = srf(isrf)%index(IND_SURFACE)

            IF (isurface.EQ.0) CYCLE

            isector  = obj(iobj)%index(IND_IS)
            isurface = isurface_list(isurface)

c            WRITE(0,*) 'SECTOR  :',iobj,isrf,isector,
c     .                 srf(isrf)%index(IND_SURFACE)

            IF (TRIM(surface(isurface)%sector).EQ.'-1') THEN
c I used to assume that a -1 mean always .TRUE. for CheckIndex,
c but this went bad when using CheckIndex with void processing,
c so I've removed this assumption...
              STOP '**** FIX -1 SPECIFICATION ****'
            ENDIF
            IF (CheckIndex(isector,surface(isurface)%sector)) CYCLE

c            STOP 'TEST'

            WRITE(0,*) 'SECTOR  :',iobj,isrf,isector,
     .                 srf(isrf)%index(IND_SURFACE)
            WRITE(0,*) '        :',isurface,isector,
     .                 srf(isrf)%index(IND_SURFACE)

            IF     (surface(isurface)%subtype.EQ.STRATUM   ) THEN
             srf(isrf)%index(IND_SURFACE) = 
     .         srf(isrf)%index(IND_SURFACE) - 1
             WRITE(0,*) 'STRATUM REMAP'
c             STOP 'TEST'
            ELSEIF (surface(isurface)%subtype.EQ.ADDITIONAL) THEN
             srf(isrf)%index(IND_SURFACE) = surface(default_surface)%num
             WRITE(0,*) 'WALL REMAP'
            ELSE
              CALL ER('ProcessTetrahedrons','Invalid surface SUBTYPE '//
     .                'when remapping the surface index',*99)
            ENDIF

            WRITE(0,*) '        :',isurface,isector,
     .                 srf(isrf)%index(IND_SURFACE)

           ENDDO
        ENDDO

c        STOP 'MADE IT HERE!'

      ENDIF

      WRITE(eirfp,*) '  NOBJ:',nobj
      WRITE(eirfp,*) '  NSRF:',nsrf
      WRITE(eirfp,*) '  NVTX:',nvtx
      WRITE(eirfp,*) 'DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: AssembleTetrahedrons
c
      SUBROUTINE AssembleTetrahedrons
      USE mod_eirene06_parameters
      USE mod_eirene06
      USE mod_eirene06_locals
      USE mod_geometry
      USE mod_options
      IMPLICIT none

      LOGICAL CheckIndex,PointInPolygon

      INTEGER itry,nseg,i1,i2,i3,itet,ivtx,iside,
     .        nlist_sec,list_sec(100),nslice,slice_i(100),
     .        nlist_sli,list_sli(100),islice,islice_index,
     .        ilist,nlist,list(100),
     .        iobj,isrf,istart,iend,fp,
     .        isurface,isurface_list(nsurface),trycycle_last
      LOGICAL default_distribution,status,hole_found
      LOGICAL :: debug  = .FALSE.
      LOGICAL :: hack   = .FALSE.
      LOGICAL :: filter = .FALSE.
      REAL*8 ang,ang1,ang2,dang(360,2),adelta,a(3,3),p(3,2),
     .       ashift
      TYPE(type_srf   ) newsrf
      TYPE(type_object) newobj

      LOGICAL, ALLOCATABLE :: trycycle(:),trycycle_global(:)
      REAL   , ALLOCATABLE :: tryxcen (:),tryycen        (:)

      WRITE(eirfp,*) 'BUILDING TETRAHEDRONS - NEW',debug,eirfp
      WRITE(0,*) 
      WRITE(0,*) ' === BUILDING TETRAHEDRONS - NEW ===',debug
      WRITE(0,*) 

      fp = 88

      IF (.TRUE.) THEN
c...    Convert legacy triangle objects to generalized geometry objects:
        ntry = ntri
        ntrysrf = 0
        ntryvtx = 0
        ALLOCATE(try   (     ntry))
        ALLOCATE(trysrf(   3*ntry))
        ALLOCATE(tryvtx(3 ,6*ntry))
        ALLOCATE(plasma(20,  ntry))  ! Needs to be consistent with TMP_PLASMA in ResolveFilament
        ALLOCATE(bfield(4 ,  ntry))
        CALL BuildNewTriangleObjects
      ENDIF

      WRITE(eirfp,*) '  NTRY:',ntry
      WRITE(eirfp,*) '  NOBJ:',nobj
      WRITE(eirfp,*) '  NSRF:',nsrf
      WRITE(eirfp,*) '  NVTX:',nvtx

c...  Start of tetrahedron objects:
      istart = nobj + 1
      
c...  Setup groups:
      ngrp = 12
      grp(3)%origin = GRP_MAGNETIC_GRID  ! *** NEED AddGroup and CopyGroup functions ***
      grp(3)%type   = GRP_TETRAHEDRON
      grp(4)%origin = GRP_VACUUM_GRID
      grp(4)%type   = GRP_TETRAHEDRON

c...  Cropping the tetrahedral grid in the poloidal plane:
      ALLOCATE(trycycle       (ntry))
      ALLOCATE(trycycle_global(ntry))
      trycycle_global = .FALSE.
      DO itet = 1, opt_eir%tet_n
        IF (opt_eir%tet_type(itet).EQ.1.0) filter = .TRUE.
      ENDDO
      IF (filter) THEN
        trycycle_global = .TRUE.
        ALLOCATE(tryxcen(ntry))
        ALLOCATE(tryycen(ntry))
        DO itry = 1, ntry
          DO i1 = 1, 3
            isrf = try(itry)%iside(i1)
            i2 = 1
            IF (isrf.LT.0) i2 = 2 ! Side orientation is switched, so use other end point
            a(1,i1) = tryvtx(1,trysrf(ABS(isrf))%ivtx(i2))
            a(2,i1) = tryvtx(2,trysrf(ABS(isrf))%ivtx(i2))
            a(3,i1) = 0.0D0
          ENDDO
          tryxcen(itry) = SNGL(SUM(a(1,1:3))) / 3.0
          tryycen(itry) = SNGL(SUM(a(2,1:3))) / 3.0
        ENDDO
        DO itet = 1, opt_eir%tet_n
          IF (opt_eir%tet_type(itet).NE.1.0) CYCLE         
          WRITE(0,*) 'ITET CROPPING:',itet
          WRITE(0,*) '    :',opt_eir%tet_type(itet)
          WRITE(0,*) '    :',opt_eir%tet_x1  (itet)
          WRITE(0,*) '    :',opt_eir%tet_y1  (itet)
          WRITE(0,*) '    :',opt_eir%tet_x2  (itet)
          WRITE(0,*) '    :',opt_eir%tet_y2  (itet)
          DO itry = 1, ntry
            IF (tryxcen(itry).GE.opt_eir%tet_x1(itet).AND.
     .          tryycen(itry).GE.opt_eir%tet_y1(itet).AND.
     .          tryxcen(itry).LE.opt_eir%tet_x2(itet).AND.
     .          tryycen(itry).LE.opt_eir%tet_y2(itet)) 
     .        trycycle_global(itry) = .FALSE.
          ENDDO
        ENDDO
        DEALLOCATE(tryxcen)
        DEALLOCATE(tryycen)
      ENDIF





c...  Toroidal distribution:

c work from the outside in:
c      -load the composite line (can I use WHERE? no I don't thing so, but how was I doing this elsewhere -- ah, from MAXVAL, or something...)
c      -need a new routine that builds a list -- BuildList -- from a character string, based on an string:
c          (1r5,2)r2 kind of thing, (1-3r2,4,1-3r2)r2 ... etc.
c      -for each sector, in turn:
c           -for each slice, in turn:
c                -determine the angular width of each slice and build up the total width of the composite
c                -also store the slice index with each slice in each section, to be used later when scanning through
c                 the slices to bluild bricks, i.e. knowing which triangles to blank out (need a tidy algorithm for this)
c

c need to reset trycycle for each slice
      default_distribution = .TRUE.

      nslice  = 0
      slice_i = 0
      dang    = 0.0D0
      adelta  = 0.0D0      

      DO itet = 1, opt_eir%tet_n
        IF (opt_eir%tet_type(itet).EQ.4.0) EXIT
      ENDDO
      IF (itet.EQ.opt_eir%tet_n+1) 
     .  CALL ER('AssembleTetrahedrons','No composite data line '//
     .          'found',*99)

      CALL GetList(nlist_sec,list_sec,opt_eir%tet_composite(itet))

      DO i1 = 1, nlist_sec
        DO itet = 1, opt_eir%tet_n
          IF (opt_eir%tet_type (itet).EQ.3.0         .AND.
     .        opt_eir%tet_index(itet).EQ.list_sec(i1)) EXIT
        ENDDO
        IF (itet.EQ.opt_eir%tet_n+1) 
     .    CALL ER('AssembleTetrahedrons','Sector data line not '//
     .            'found',*99)
        CALL GetList(nlist_sli,list_sli,opt_eir%tet_sec_list(itet))
      
        DO i2 = 1, nlist_sli
          DO itet = 1, opt_eir%tet_n
            IF (opt_eir%tet_type (itet).EQ.2.0         .AND.
     .          opt_eir%tet_index(itet).EQ.list_sli(i2)) EXIT
          ENDDO
          IF (itet.EQ.opt_eir%tet_n+1) 
     .      CALL ER('AssembleTetrahedrons','Slice data line not '//
     .              'found',*99)

          WRITE(0,*) 'processing:',itet,opt_eir%tet_index(itet)

          nslice = nslice + 1
          slice_i(nslice) = itet  
          SELECTCASE(opt_eir%tet_mode(itet))
            CASE (1)
              dang(nslice,1) = 0.0D0
              dang(nslice,2) = opt_eir%tet_param1(itet) 
              IF (dang(nslice,1).GT.dang(nslice,2))
     .          CALL ER('AssembleTetrahedrons','Bad angle range',*99)
              adelta = adelta + dang(nslice,1) 
            CASE DEFAULT
              CALL ER('AssembleTetrahedrons','Unknown mode',*99)
          ENDSELECT

          default_distribution = .FALSE.
        ENDDO ! Slices
      ENDDO ! Sectors

      WRITE(0,*) 'SLICE INDEX: '
      DO i1 = 1, nslice
        WRITE(0,*) i1,slice_i(i1),dang(i1,1:2)
      ENDDO



      IF (default_distribution) THEN 
        ang = 360.0D0 / DBLE(ntorseg) * torfrac
        nslice = ntorseg
        dang(:,1) = 0.0D0
        dang(:,2) = ang
        DO islice = 1, nslice
          adelta = adelta + dang(islice,2) - dang(islice,1) 
        ENDDO
c        dang(1,1) = 0.0D0
c        dang(1,2) = dang(1,1) + ang
c        DO islice = 2, nslice
c          dang(islice,1) = dang(islice-1,2)
c          dang(islice,2) = dang(islice  ,1) + ang
c        ENDDO
c        adelta = dang(nslice,2) - dang(1,1)
c        dang = dang - 0.5D0 * adelta      
      ENDIF

      ashift = 0.0D0
      DO islice = 1, nslice
        WRITE(0,'(A,I6,2F12.6)') 
     .    'ANGLES:',islice,dang(islice,1:2)+ashift
        WRITE(eirfp,'(A,I6,2F12.6)') 
     .    'ANGLES:',islice,dang(islice,1:2)
        ashift = ashift + dang(islice,2)
      ENDDO



      dang = dang * D_DEGRAD
      adelta = adelta * D_DEGRAD


c...  
      ashift = 0.0D0
      trycycle_last = -1
      DO islice = 1, nslice

        islice_index = slice_i(islice)

        WRITE(0 ,*) '===== tet slice ===>',islice_index,
     .              opt_eir%tet_index(islice_index)
        WRITE(88,*) '===== tet slice ===>',islice_index

        IF (trycycle_last.NE.islice_index) THEN
          trycycle = .FALSE.

          CALL GetList(nlist,list,opt_eir%tet_del_hole(islice_index))
          WRITE(0,*) '----->list',nlist,' '//
     .       TRIM(opt_eir%tet_del_hole(islice_index))

          hole_found = .TRUE.
          IF (nlist.GT.0) hole_found = .FALSE.

          DO ilist = 1, nlist
c            WRITE(0,*) '      ----> go ',list(ilist)
            DO isurface = 1, nsurface            
              IF (surface(isurface)%type    .EQ.HOLE_IN_GRID.AND.
     .            surface(isurface)%index(2).EQ.list(ilist)) EXIT
            ENDDO
            IF (isurface.EQ.nsurface+1) CYCLE

            hole_found = .TRUE.
            
            DO i1 = 1, ntry
              DO i2 = 1, 3
                isrf = try(i1)%iside(i2)
                i3 = 1
                IF (isrf.LT.0) i3 = 2 ! Side orientation is switched, so use other end point
                p(i2,1) = tryvtx(1,trysrf(ABS(isrf))%ivtx(i3))
                p(i2,2) = tryvtx(2,trysrf(ABS(isrf))%ivtx(i3))
              ENDDO
              IF (PointInPolygon(surface(isurface)%v(1,1),
     .                           surface(isurface)%v(2,1),
     .                           3,p(1,1))) THEN
c                WRITE(0,*) 'FOUND!',i1
c                WRITE(0,*) '    t1',p(1,1:2)
c                WRITE(0,*) '    t2',p(2,1:2)
c                WRITE(0,*) '    t3',p(3,1:2)
c                WRITE(0,*) '      ',surface(isurface)%v(1:2,1)

                EXIT

              ENDIF 
            ENDDO
            IF (i1.EQ.itry+1) 
     .        CALL ER('AssembleTetrahedrons','Hole triangle '//
     .                'not found',*99)

            trycycle(i1) = .TRUE.

            status = .TRUE.
            DO WHILE (status)
c              WRITE(0,*) '        ----> pass '
              status = .FALSE.
              DO i1 = 1, ntry
                IF (.NOT.trycycle(i1)) CYCLE
                DO i2 = 1, 3
                  isrf = ABS(try(i1)%iside(i2))
                  IF (trysrf(isrf)%index(IND_SURFACE).NE.0) CYCLE
                  i3 = try(i1)%omap(i2)
                  IF (.NOT.trycycle(i3)) THEN
                    trycycle(i3) = .TRUE.
                    status = .TRUE.
c                    WRITE(0,*) '  blocking',i1,
c     .           tryvtx(1,trysrf(ABS(isrf))%ivtx(1)),
c     .           tryvtx(2,trysrf(ABS(isrf))%ivtx(1))
                  ENDIF
                ENDDO 
              ENDDO
            ENDDO

          ENDDO

          IF (.NOT.hole_found)
     .      CALL ER('AssembleTetrahedrons','Hole not found',*99)

          trycycle = trycycle.OR.trycycle_global
        ENDIF



        ang1 = dang(islice,1) + ashift
        ang2 = dang(islice,2) + ashift

        WRITE(eirfp,'(A,I6,2F12.6)') '   TOROIDALIZING:',islice,
     .    ang1/D_DEGRAD,ang2/D_DEGRAD
     .    

c left off
c check the stratum labels in .eirdat and figure out why they don't match whatś in the in put files
c compare a case with a transparent target and one without
c ...there was something else but i cant remember what...

        CALL BuildBricks(islice,opt_eir%tet_index(islice_index),
     .                   ang1,ang2,trycycle)

        WRITE(eirfp,*) '  NTRY:',ntry
        WRITE(eirfp,*) '  NOBJ:',nobj
        WRITE(eirfp,*) '  NSRF:',nsrf
        WRITE(eirfp,*) '  NVTX:',nvtx

        ashift = ashift + dang(islice,2)

        trycycle_last = islice_index






c        isurface_list = 0
c        DO i1 = 1, nsurface
c          IF (surface(i1)%type.NE.NON_DEFAULT_STANDARD) CYCLE
c          isurface_list(surface(i1)%num) = i1
c        ENDDO
c        DO iobj = 1, nobj
c          DO iside = 1, obj(iobj)%nside
c            isrf     = ABS(obj(iobj)%iside(iside))
c            isurface = srf(isrf)%index(IND_SURFACE)
c            IF (isurface.EQ.0) CYCLE
c            isurface = isurface_list(isurface)
c
c            WRITE(0 ,*) 'slicing check',islice,
c     .                   opt_eir%tet_index(islice_index),
c     .                   obj(iobj)%index(IND_ISI),isurface
c
c            WRITE(88,*) 'slicing check',islice,
c     .                   opt_eir%tet_index(islice_index),
c     .                   obj(iobj)%index(IND_ISI),isurface
c          ENDDO
c        ENDDO

      ENDDO

c...  Done with these, clear some memory:
      DEALLOCATE(trycycle)
      DEALLOCATE(try)
      DEALLOCATE(trysrf)
      DEALLOCATE(tryvtx)
 


C     IN RAY/OUT, I DON'T YET HAVE A SITUATION WHERE ORIENTATION IS A BIG ISSUE, SO FOR EIRENE
c     I'LL JUST DO THE ORDERING ON THE FLY, CHECKING IF THE POINTS ARE CO- OR COUNTER CLOCKWISE
c     WITH RESPECT TO THE CENTER OF THE TETRAHEDRON, BUT NEED TO COME UP WITH A TEST FIRST, 
c     PERHAPS THE NORMAL WITH RESPECT TO THE CENTER VECTOR? 

      iend = nobj

      WRITE(eirfp,*) '  NOBJ:',nobj
      WRITE(eirfp,*) '  NSRF:',nsrf
      WRITE(eirfp,*) '  NVTX:',nvtx

c...  Build connection map:
      CALL BuildConnectionMap(istart,iend)

c.... For Eirene, need to make sure the tetrahedrons are arranged according to
c     convention:
      CALL FixTetrahedrons(istart,iend)




c        isurface_list = 0
c        DO i1 = 1, nsurface
c          IF (surface(i1)%type.NE.NON_DEFAULT_STANDARD) CYCLE
c          isurface_list(surface(i1)%num) = i1
c        ENDDO
c        DO iobj = 1, nobj
c          DO iside = 1, obj(iobj)%nside
c            isrf     = ABS(obj(iobj)%iside(iside))
c            isurface = srf(isrf)%index(IND_SURFACE)
c            IF (isurface.EQ.0) CYCLE
c            isurface = isurface_list(isurface)
c            WRITE(0 ,*) 'slicing check 2',
c     .                   obj(iobj)%index(IND_ISI),isurface
c            WRITE(88,*) 'slicing check 2',
c     .                   obj(iobj)%index(IND_ISI),isurface
c          ENDDO
c        ENDDO





c...  Plasma association:
c      CALL CheckTetrahedronStructure  ! Add this here..?

c...  Add slices here, if requested:

c...  Impose filament structures:
      IF (opt_fil%opt.NE.0) THEN
        WRITE(eirfp,*) '  NOBJ:',nobj
        WRITE(eirfp,*) '  NSRF:',nsrf
        WRITE(eirfp,*) '  NVTX:',nvtx
        CALL ResolveFilament(-1)
        CALL AssignFilamentPlasma   
        iend = nobj
        CALL BuildConnectionMap(istart,iend)
        CALL FixTetrahedrons(istart,iend)
      ENDIF




c...  Manual refinement:
      DO itet = 1, opt_eir%tet_n
        IF (opt_eir%tet_type(itet).NE.5.0) CYCLE
        CALL ResolveFilament(itet)
        iend = nobj
        CALL BuildConnectionMap(istart,iend)
        CALL FixTetrahedrons   (istart,iend)
      ENDDO

c      isurface_list = 0
c      DO i1 = 1, nsurface
c        IF (surface(i1)%type.NE.NON_DEFAULT_STANDARD) CYCLE
c        isurface_list(surface(i1)%num) = i1
c      ENDDO
c      DO iobj = 1, nobj
c        DO iside = 1, obj(iobj)%nside
c          isrf     = ABS(obj(iobj)%iside(iside))
c          isurface = srf(isrf)%index(IND_SURFACE)
c          IF (isurface.EQ.0) CYCLE
c          isurface = isurface_list(isurface)
c          WRITE(0 ,*) 'slicing check 3',
c     .                 obj(iobj)%index(IND_ISI),isurface
c          WRITE(88,*) 'slicing check 3',
c     .                 obj(iobj)%index(IND_ISI),isurface
c        ENDDO
c      ENDDO



      CALL CheckTetrahedronStructure






c      isurface_list = 0
c      DO i1 = 1, nsurface
c        IF (surface(i1)%type.NE.NON_DEFAULT_STANDARD) CYCLE
c        isurface_list(surface(i1)%num) = i1
c      ENDDO
c      DO iobj = 1, nobj
c        DO iside = 1, obj(iobj)%nside
c          isrf     = ABS(obj(iobj)%iside(iside))
c          isurface = srf(isrf)%index(IND_SURFACE)
c          IF (isurface.EQ.0) CYCLE
c          isurface = isurface_list(isurface)
c          WRITE(0 ,*) 'slicing check 4',
c     .                 obj(iobj)%index(IND_ISI),isurface
c          WRITE(88,*) 'slicing check 4',
c     .                 obj(iobj)%index(IND_ISI),isurface
c        ENDDO
c      ENDDO



      IF (.TRUE.) THEN
c...    Make a linked list for the non-default standard surfaces in the surface 
c       array (which also includes wall surfaces, which are not of interest here):
        isurface_list = 0
        DO i1 = 1, nsurface
          IF (surface(i1)%type.NE.NON_DEFAULT_STANDARD) CYCLE
          isurface_list(surface(i1)%num) = i1
        ENDDO

c        WRITE(fp,*) 'LIST:'
c        WRITE(fp,*) isurface_list(1:nsurface)

        DO iobj = 1, nobj
c...      Collect connection map information:
          DO iside = 1, obj(iobj)%nside
            isrf   = ABS(obj(iobj)%iside(iside))

            isurface = srf(isrf)%index(IND_SURFACE)

            IF (isurface.EQ.0) CYCLE

            isurface = isurface_list(isurface)
            islice   = obj(iobj)%index(IND_ISI)


c            WRITE(fp,*) 'SECTOR 1:',iobj,isrf,islice,isurface,
c     .                 ' '//TRIM(surface(isurface)%sector)

            IF (TRIM(surface(isurface)%sector).EQ.'-1') THEN
c I used to assume that a -1 mean always .TRUE. for CheckIndex,
c but this went bad when using CheckIndex with void processing,
c so I've removed this assumption...
              STOP '**** FIX -1 SPECIFICATION ****'
            ENDIF



            IF (surface(isurface)%sector(1:3).NE.'all') THEN
c              WRITE(0 ,*) 'SECTOR 2:',iobj,isrf,obj(iobj)%index(IND_IS),
c     .                    islice,' '//TRIM(surface(isurface)%sector)
              WRITE(88,*) 'SECTOR 2:',iobj,isrf,obj(iobj)%index(IND_IS),
     .                    islice,' '//TRIM(surface(isurface)%sector)
            ENDIF

            IF (islice.EQ.2) THEN
c              WRITE(0 ,*) 'slicing',islice,isurface
              WRITE(88,*) 'slicing',islice,isurface
            ENDIF

c...        This check seems backwards, and so it can hurt the brain, but
c           it works because the surface properties over-ride is applied to
c           all slices by default, so the objective here is to identify
c           slices that are not included in the over-ride and to point them
c           to some default surface:
            IF (CheckIndex(islice,surface(isurface)%sector)) CYCLE



            WRITE(fp,*) 'SECTOR 2:',iobj,isrf,obj(iobj)%index(IND_IS),
     .                  islice,
     .                  ' '//TRIM(surface(isurface)%sector)
            WRITE(fp,*) '        :',isurface,islice,
     .                  srf(isrf)%index(IND_SURFACE)

            IF     (surface(isurface)%subtype.EQ.STRATUM   ) THEN
              srf(isrf)%index(IND_SURFACE) = 
     .                                  srf(isrf)%index(IND_SURFACE) - 1
             WRITE(fp,*) 'STRATUM REMAP'
c             STOP 'TEST'
            ELSEIF (surface(isurface)%subtype.EQ.ADDITIONAL) THEN
             srf(isrf)%index(IND_SURFACE) = surface(default_surface)%num
             WRITE(fp,*) 'WALL REMAP'
            ELSE
              CALL ER('ProcessTetrahedrons','Invalid surface SUBTYPE '//
     .                'when remapping the surface index',*99)
            ENDIF

            WRITE(fp,*) '        :',isurface,islice,
     .                  srf(isrf)%index(IND_SURFACE)

           ENDDO
        ENDDO

      ENDIF


c     Looking for the special tetrahedron catch all surface of desperation...
      IF (surface(nsurface)%type    .NE.NON_DEFAULT_STANDARD.OR.
     .    surface(nsurface)%subtype .NE.ADDITIONAL          .OR.
     .    surface(nsurface)%index(1).NE.-1) 
     .  CALL ER('AssembleTetrahedrons','Dump surface not found',*99)

      DO iobj = istart, iend
        DO iside = 1, obj(iobj)%nside
          isrf = ABS(obj(iobj)%iside(iside))
          IF (iside.EQ.1.AND.
     .        obj(iobj)%omap(iside)       .EQ.0.AND.
     .        srf(isrf)%index(IND_SURFACE).EQ.0) THEN

            srf(isrf)%index(IND_SURFACE) = surface(nsurface)%num

          ENDIF

          IF (iside.EQ.1.AND.
     .        obj(iobj)%omap(iside)       .EQ.0.AND.
     .        srf(isrf)%index(IND_SURFACE).NE.0) THEN

            DO i1 = 1, nsurface
              IF (surface(i1)%type.EQ.NON_DEFAULT_STANDARD.AND.
     .            surface(i1)%num .EQ.srf(isrf)%index(IND_SURFACE))
     .          EXIT
            ENDDO
            IF (i1.EQ.nsurface+1) THEN
              WRITE(0,*) '  PROBLEM INDEX=',srf(isrf)%index(IND_SURFACE)
              DO i1 = 1, nsurface
                IF (surface(i1)%type.NE.NON_DEFAULT_STANDARD) CYCLE
                WRITE(0,*) '  INDICES=',i1,surface(i1)%index(5),
     .                       surface(i1)%num
              ENDDO
              CALL ER('AssembleTetrahedrons','Surface not found',*99)
            ENDIF

            IF (surface(i1)%iliin.LE.0) THEN
c              WRITE(0,*) 'cleaning',obj(iobj)%index(IND_IS ),
c     .                              obj(iobj)%index(IND_ISI)
              srf(isrf)%index(IND_SURFACE) = surface(nsurface)%num
            ENDIF

          ENDIF

          IF (iside.NE.1.AND.
     .        obj(iobj)%omap(iside)       .EQ.0.AND.
     .        srf(isrf)%index(IND_SURFACE).LE.0) THEN
            WRITE(0,*) 'SOMETHING IS WRONG'
            STOP
          ENDIF

        ENDDO
      ENDDO


      WRITE(eirfp,*) '  NOBJ:',nobj
      WRITE(eirfp,*) '  NSRF:',nsrf
      WRITE(eirfp,*) '  NVTX:',nvtx
      WRITE(eirfp,*) 'DONE'



      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: BuildBricks
c
      SUBROUTINE BuildBricks(islice,islice_index,ang1,ang2,trycycle)
      USE mod_eirene06_parameters
      USE mod_eirene06
      USE mod_eirene06_locals
      USE mod_geometry
      USE mod_options
      IMPLICIT none

      INTEGER, INTENT(IN)  :: islice,islice_index
      LOGICAL, INTENT(IN)  :: trycycle(*)
      REAL*8 , INTENT(IN)  :: ang1,ang2

      INTEGER itry,itry1,i1,i2,isid,itet,ivtx,iobj,iside,isrf,ishift,
     .        idum1
      LOGICAL :: debug  = .FALSE.
      LOGICAL :: hack   = .FALSE.
      LOGICAL :: trycheck(ntry)
      REAL*8 a(3,3),b(3,7),c(3,4)
      TYPE(type_srf   ) newsrf
      TYPE(type_object) newobj

      REAL*8     DTOL
      PARAMETER (DTOL=1.0D-07)

      INTEGER s(4,14)
      DATA    s / 3, 2, 1, 7,   
     .            1, 2, 5, 7,   5, 4, 1, 7,   
     .            2, 3, 6, 7,   6, 5, 2, 7,   
     .            3, 1, 4, 7,   4, 6, 3, 7,   
     .            4, 5, 6, 7,
c                Reversed:
     .            1, 2, 4, 7,   4, 2, 5, 7,
     .            2, 3, 5, 7,   5, 3, 6, 7,
     .            3, 1, 6, 7,   6, 1, 4, 7 /

      INTEGER t(3,4)
      DATA    t / 1, 2, 3,  1, 4, 2,  2, 4, 3,  3, 4, 1 /


      WRITE(eirfp,*) 'BUILDING TETRAHEDRONS - NEW',debug,eirfp

      trycheck = .FALSE.

      DO itry = 1, ntry  ! ntry   

        IF (trycycle(itry)) CYCLE

c...    Assemble brick vertices:
        DO i1 = 1, 3
          isrf = try(itry)%iside(i1)
          i2 = 1
          IF (isrf.LT.0) i2 = 2 ! Side orientation is switched, so use other end point
          a(1,i1) = tryvtx(1,trysrf(ABS(isrf))%ivtx(i2))
          a(2,i1) = tryvtx(2,trysrf(ABS(isrf))%ivtx(i2))
          a(3,i1) = 0.0D0
        ENDDO

c...    Check orientation...?
        IF (debug) THEN
          WRITE(eirfp,*) 'A:',a(1:2,1)
          WRITE(eirfp,*) 'A:',a(1:2,2)
          WRITE(eirfp,*) 'A:',a(1:2,3)
        ENDIF

c...    Expand toroidally:

c        isector = 1
        IF (debug) THEN
          WRITE(eirfp,*) 'DANG:',ang1/D_DEGRAD,ang2/D_DEGRAD
        ENDIF
        b(1,1:3) = a(1,1:3) * DCOS(ang1)  !DCOS(dang(isector,1))
        b(2,1:3) = a(2,1:3)
        b(3,1:3) = a(1,1:3) * DSIN(ang1)  !DSIN(dang(isector,1))
        b(1,4:6) = a(1,1:3) * DCOS(ang2)  !DCOS(dang(isector,2))
        b(2,4:6) = a(2,1:3)
        b(3,4:6) = a(1,1:3) * DSIN(ang2)  !DSIN(dang(isector,2))

c       Center of brick:
        b(1:3,7) = 0.0D0
        DO i1 = 1, 6
          b(1:3,7) = b(1:3,7) + b(1:3,i1) / 6.0D0
        ENDDO

c...    Check orientation...?
        IF (debug) THEN
          WRITE(eirfp,*) 'B:',b(1:3,1)
          WRITE(eirfp,*) 'B:',b(1:3,2)
          WRITE(eirfp,*) 'B:',b(1:3,3)
          WRITE(eirfp,*) 'B:',b(1:3,4)
          WRITE(eirfp,*) 'B:',b(1:3,5)
          WRITE(eirfp,*) 'B:',b(1:3,6)
          WRITE(eirfp,*) 'B:',b(1:3,7)
        ENDIF

c...    Determine group assignment (ad hoc at the moment):

c...    Make tetrahedrons:
        DO itet = 1, 8  ! 8 tetrahedrons for each triangle
          newobj = try(itry)
          newobj%segment(1) = itet
          newobj%group = try(itry)%group + 2 
          newobj%phi = SNGL(0.5D0*(ang1 + ang2) / D_DEGRAD)
c          newobj%phi = SNGL(0.5D0*(dang(1,1) + dang(1,2)) / D_DEGRAD)
          newobj%nside = 4
          newobj%index(IND_IS ) = islice ! isector ! 1
          newobj%index(IND_ISI) = islice_index
          newobj%index(IND_IK) = try(itry)%index(IND_IK)
          newobj%index(IND_IR) = try(itry)%index(IND_IR)

          ishift = 0
          IF (itet.GE.2.AND.itet.LE.7) THEN
            iside = INT(REAL(itet) / 2.0)
            itry1 = try(itry)%omap(iside)
            IF (itry1.NE.0) THEN
              IF (trycheck(itry1)) ishift = 7
            ENDIF
          ENDIF

          c(1:3,1) = b(1:3,s(1,itet+ishift))
          c(1:3,2) = b(1:3,s(2,itet+ishift))
          c(1:3,3) = b(1:3,s(3,itet+ishift))
          c(1:3,4) = b(1:3,s(4,itet+ishift))

          IF (debug) THEN
            WRITE(eirfp,*)
            WRITE(eirfp,*) itet,ishift
            WRITE(eirfp,*) s(1:4,itet+ishift)
            WRITE(eirfp,*) 'C:',c(1:3,1)
            WRITE(eirfp,*) 'C:',c(1:3,2)
            WRITE(eirfp,*) 'C:',c(1:3,3)
            WRITE(eirfp,*) 'C:',c(1:3,4)
          ENDIF

c...      Assign index mapping:
          DO isid = 1, newobj%nside
c           Vertices:
            IF (isid.EQ.1.AND.itet.GE.2.AND.itet.LE.7) THEN
              i1 = INT(REAL(itet-1)/2.0+0.51)                
              IF (debug) WRITE(eirfp,*) 'IIII:',itet,i1
              newsrf%index = trysrf(ABS(try(itry)%iside(i1)))%index
          
              IF (debug) THEN
                WRITE(eirfp,*) ' :',
     .            trysrf(ABS(try(itry)%iside(i1)))%index(1:3)
              ENDIF
            ELSE
              newsrf%index = 0
            ENDIF
            newsrf%type = SPR_PLANAR_POLYGON
            newsrf%obj  = nobj + 1
            newsrf%side = isid
            newsrf%nvtx = 3
            DO ivtx = 1, 3
c              IF (DABS(c(1,t(ivtx,isid))).LT.1.0D-07) THEN
c                WRITE(fp,*) 'TETRAHEDRONS: SUSPICIOUS X-VAL A',nobj+1
c                STOP 
c              ENDIF
              newsrf%ivtx(ivtx) = AddVertex(c(1,t(ivtx,isid)))
            ENDDO
            
            newobj%iside(isid) = AddSurface(newsrf)

            IF (debug) THEN
              WRITE(eirfp,*) 'I.:',isid
              WRITE(eirfp,*) 'IN:',newsrf%index(1:3)
              WRITE(eirfp,*) 'I0:',t(1:3,isid)
              WRITE(eirfp,*) 'I1:',newsrf%ivtx(1:3)
              WRITE(eirfp,*) 'I1:',srf(ABS(newobj%iside(isid)))%
     .                                     ivtx(1:3)

              WRITE(eirfp,*) 'V:',vtx(1:2,newsrf%ivtx(1))
              WRITE(eirfp,*) 'V:',vtx(1:2,newsrf%ivtx(2))
              WRITE(eirfp,*) 'V:',vtx(1:2,newsrf%ivtx(3))
            ENDIF

          ENDDO
c...      Center of tetrahedron -- needs to be done properly:
          newobj%x = SNGL(0.25D0 * SUM(c(1,1:4)))
          newobj%y = SNGL(0.25D0 * SUM(c(2,1:4)))
          newobj%z = SNGL(0.25D0 * SUM(c(3,1:4)))

c          IF (itry.EQ.1) THEN
c            WRITE(fp,*) 'Cx:',c(1,1:4)
c            WRITE(fp,*) 'Cy:',c(2,1:4)
c            WRITE(fp,*) 'Cz:',c(3,1:4)
c            WRITE(fp,*) 'Cen:',newobj%x,newobj%y,newobj%z
c          ENDIF

c...      Add object:
          idum1 = AddObject(newobj)

        ENDDO  ! Tetrahedron loop

        trycheck(itry) = .TRUE.
         
      ENDDO  ! Triangle loop


      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: ProcessTriangles_06
c
      SUBROUTINE ProcessTriangles_06(mode)
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER mode


      LOGICAL PointOnLine
      REAL*8  MaxTriangleAngle,TriangleSideLength

      REAL       TOL        ,DTOL
c      PARAMETER (TOL=1.0E-05,DTOL=1.0D-07)
      PARAMETER (TOL=1.0E-06,DTOL=1.0D-06)

      INTEGER i1,i2,i3,i4,v1,v2,v3,v4,knot,ring,side,target,
     .        xupdate(10),yupdate(10),ix,iy,iscan,problem_triangle
      LOGICAL test,output,malformed,dummy_test, ! surface_assigned,
     .        wall_assignment,grid_assignment,warning_given
      REAL    xmin,xmax,ymin,ymax,xval,yval
      REAL*8  x(0:2),y(0:2),s,t,dist1,dist2,dist3,dist4

      INTEGER, ALLOCATABLE :: xregion(:),yregion(:),nregion(:,:),
     .                        iregion(:,:,:)

      DATA (xupdate(i1),i1=1,10) /0, -1,  0,  1, -1, 1, -1, 0, 1, 0/ , 
     .     (yupdate(i1),i1=1,10) /0, -1, -1, -1,  0, 0,  1, 1, 1, 0/

      warning_given = .FALSE.

      eir_pass = eir_pass + 1
  
      WRITE(eirfp,*) 'PROCESSING TRIANGLES'  

      output = .FALSE.

      problem_triangle = 893

      IF (mode.EQ.-1) GOTO 10

      WRITE(eirfp,*) '  REMOVING DUPLICATE VERTICIES'

c...  Check if there are any malformed triangles:
      malformed = .FALSE.
      DO i1 = ntri, 1, -1
        DO v1 = 1, 3   
          v2 = v1 + 1
          IF (v1.EQ.3) v2 = 1
          IF (tri(i1)%ver(v1).EQ.tri(i1)%ver(v2)) THEN
            WRITE(eirfp,*) 'MALFORMED TRIANGLE DETECTED',i1
            WRITE(eirfp,*) 'IK,IR=',tri(i1)%index(1:2)
            malformed = .TRUE.
          ENDIF
        ENDDO
      ENDDO
      IF (malformed) THEN
c       This was being triggered with the ITER grid iterm.carre.105 because some cells
c       near IRWALL were being incorrectly registered as inside the current
c       focus cell in ProcessFluidGrid, thanks to a lax DTOL in PointOnLine.  PointOnLine is a
c       problem routine but I can't figure out a better solution as long as DIVIMP
c       uses single precision REAL with RVERTP, etc., i.e. I need to use some DTOL threshold.
c       Best if everything were REAL*8. - SL, 23.07.09 
        CALL WriteEireneTriangles
        CALL SaveTriangles_06
        CALL DumpGrid('MALFORMED FLUID GRID TRIANGLES FOUND')
      ENDIF

c...  Eliminate duplicate verticies:    ! SPEED:? SORT VERTICIES INTO REGIONS AND ONLY SCAN OVER NEIGHBOUR REGIONS? 
      DO i1 = 1, nver
        DO i2 = i1+1, nver
          IF (ver(i1,1).NE.-999.0D0.AND.
     .        DABS(ver(i1,1)-ver(i2,1)).LT.DTOL.AND.
     .        DABS(ver(i1,2)-ver(i2,2)).LT.DTOL) THEN
c            WRITE(eirfp,*) '---BYE---'
c            WRITE(eirfp,*) i1,i2,ver(i1,1),ver(i1,2)
c            WRITE(eirfp,*) i1,i2,ver(i2,1),ver(i2,2)
            ver(i2,1) = -999.0D0
            ver(i2,2) = -999.0D0
c          IF (ver(i1,1).NE.-999.0.AND.
c     .        ABS(ver(i1,1)-ver(i2,1)).LT.TOL.AND.
c     .        ABS(ver(i1,2)-ver(i2,2)).LT.TOL) THEN
c            ver(i2,1) = -999.0
c            ver(i2,2) = -999.0
c...        Search all triangles for the vertex to be removed and update:
            DO i3 = 1, ntri
              DO v1 = 1, 3
                IF (tri(i3)%ver(v1).EQ.i2) tri(i3)%ver(v1) = i1
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
c...  Remove vertices tagged for deletion:
      DO i1 = nver, 1, -1
        IF (ver(i1,1).EQ.-999.0D0) THEN
c        IF (ver(i1,1).EQ.-999.0) THEN
          DO i2 = i1, nver-1
            ver(i2,1) = ver(i2+1,1)
            ver(i2,2) = ver(i2+1,2)
          ENDDO
          DO i2 = 1, ntri
            DO v1 = 1, 3
              IF (tri(i2)%ver(v1).GE.i1) 
     .          tri(i2)%ver(v1) = tri(i2)%ver(v1) - 1
            ENDDO
          ENDDO
          nver = nver - 1
        ENDIF
      ENDDO

c...  Check if there are any malformed triangles (again):
      malformed = .FALSE.
      DO i1 = ntri, 1, -1
        DO v1 = 1, 3   
          v2 = v1 + 1
          IF (v1.EQ.3) v2 = 1
          IF (tri(i1)%ver(v1).EQ.tri(i1)%ver(v2)) THEN
            WRITE(eirfp,*) 'MALFORMED TRIANGLE DETECTED',i1
            WRITE(eirfp,*) 'IK,IR=',tri(i1)%index(1:2)
            malformed = .TRUE.
          ENDIF
        ENDDO
      ENDDO
      IF (malformed) THEN
        CALL WriteEireneTriangles
        CALL SaveTriangles_06
        CALL DumpGrid('MALFORMED FLUID GRID TRIANGLES FOUND #2')
      ENDIF

cc...  Check if 2 close lying points got merged by accident, and if 
cc     yes, then remove the triangle:
c      DO i1 = ntri, 1, -1
c        DO v1 = 1, 3   
c          v2 = v1 + 1
c          IF (v1.EQ.3) v2 = 1
c          IF (tri(i1)%ver(v1).EQ.tri(i1)%ver(v2)) THEN
c            WRITE(eirfp,*) 'KILLING TRIANLGE-FIX!:',i1
c            DO i2 = i1, ntri-1
c              tri(i2) = tri(i2+1)
c            ENDDO
c            ntri = ntri - 1
c            EXIT
c          ENDIF
c        ENDDO
c      ENDDO

 10   CONTINUE

c...  (Lame) integrity checks
      DO i1 = 1, ntri
        IF (MaxTriangleAngle(i1).GT.179.0D0) 
     .    WRITE(eirfp,'(1X,A,I6,F10.2)') 
     .      'WARNING: Large interior angle detected for '//
     .      'triangle  ',i1,MaxTriangleAngle(i1)
        DO i2 = 1, 3
          IF (TriangleSideLength(i1,i2).LT.1.0D-06) 
     .    WRITE(eirfp,'(1X,A,2I6,1P,2E10.2,0P)') 
     .      'WARNING: Very short side detected for '//
     .      'triangle  ',i1,i2,TriangleSideLength(i1,i2),
     .      MaxTriangleAngle(i1)
        ENDDO
      ENDDO

      WRITE(eirfp,*) '  ASSIGNING REGIONS' 
      ALLOCATE(xregion(ntri)) ! Move to triangle construct 
      ALLOCATE(yregion(ntri))
      ALLOCATE(nregion(10,10)) 
      ALLOCATE(iregion(10,10,ntri))  ! I used to divide NTRI by 2 or 3 to reduce the array size, but this causes
                                     ! problems if the grid is very dense in a particular region of interest 

c...  Find horizontal and vertical extents of the grid:
      xmin =  1.0E+31
      xmax = -1.0E+31
      ymin =  1.0E+31
      ymax = -1.0E+31
      DO i1 = 1, ntri      
        DO v1 = 1, 3
          xmin = MIN(xmin,SNGL(ver(tri(i1)%ver(v1),1)))
          xmax = MAX(xmax,SNGL(ver(tri(i1)%ver(v1),1)))
          ymin = MIN(ymin,SNGL(ver(tri(i1)%ver(v1),2)))
          ymax = MAX(ymax,SNGL(ver(tri(i1)%ver(v1),2)))
        ENDDO
      ENDDO
      xmin = xmin - 0.01
      xmax = xmax + 0.01
      ymin = ymin - 0.01
      ymax = ymax + 0.01
c...  Setup bins:
      v1 = 1
      DO i1 = 1, ntri      
        xval = SNGL(ver(tri(i1)%ver(v1),1))  ! ***BETTER TO USE TRIANGLE CENTERS, BUT NOT STORED YET... 
        yval = SNGL(ver(tri(i1)%ver(v1),2))
        xregion(i1) = INT((xval - xmin) / (xmax - xmin) * 10.0) + 1
        yregion(i1) = INT((yval - ymin) / (ymax - ymin) * 10.0) + 1
      ENDDO
c...  Build list of regions:
      WRITE(eirfp,*) '  BUILDING REGION LISTS' 
      nregion = 0
      iregion = 0
      DO i1 =  1, ntri
        ix = xregion(i1)
        iy = yregion(i1)
        nregion(ix,iy)                = nregion(ix,iy) + 1
        iregion(ix,iy,nregion(ix,iy)) = i1
      ENDDO
      DO i1 = 10, 1, -1
        WRITE(eirfp,'(2X,10I5)') (nregion(i2,i1),i2=1,10)
      ENDDO

      WRITE(eirfp,*) '  BUILDING CONNECTION MAP',ntri

c...  Build connection map:
      DO i1 = 1, ntri
        tri(i1)%map(1:3) = 0
        tri(i1)%sid(1:3) = 0
        DO v1 = 1, 3

          v2 = v1 + 1
          IF (v1.EQ.3) v2 = 1

c...      Find matching cell based on vertex coordinates:
          ix = xregion(i1)
          iy = yregion(i1)
          iscan = 0
          i3 = 0
          DO WHILE (iscan.LE.10)

c...        Advance the loop index 'i2' (sorry that this is convoluted, but it needs 
c           to be fast 'cause some grids are big):
            i3 = i3 + 1
            IF     (iscan.EQ.0) THEN
c...          Check if a map is already assigned, and if yes, look there 
c             first (iscan=0):
              IF (tri(i1)%map(v1).GT.0) THEN
                i2 = tri(i1)%map(v1)
              ELSE
                i3 = 0
                iscan = 1
                CYCLE
              ENDIF
            ELSEIF (iscan.EQ.10) THEN
c...          Map not found so search the whole triangle mesh to be sure that 
c             there isn't one -- a large triangle can cause problems for 0<iscan<10 (iscan=10):
c              IF (i1.EQ.423) WRITE(0,*) ' WHF AS:',i2,i3,iscan
              IF (i2.LT.ntri) THEN
                i2 = i3
c                IF (i1.EQ.423) WRITE(0,*) ' WHF:',i2,iscan
              ELSE
c...            All triangles have been searched but no mapping found, trigger exit condition:
                iscan = 999
                CYCLE
              ENDIF
            ELSE
c...          First, search the region that the triangle is in (iscan=1), and 
c             then the neighbouring regions:
              IF (i3.LE.nregion(ix,iy)) THEN
                i2 = iregion(ix,iy,i3)
              ELSE 
c...            Change search region, making sure the new region is valid:
                i3 = 0
                i2 = 0
c                IF (i1.EQ.423) WRITE(0,*) ' WHF ??:',i3,iscan
                ix = 0
                DO WHILE (ix.LT.1.OR.ix.GT.10.OR.
     .                    iy.LT.1.OR.iy.GT.10)
                  iscan = iscan + 1
                  ix = xregion(i1) + xupdate(iscan)
                  iy = yregion(i1) + yupdate(iscan)
                ENDDO        
                CYCLE
              ENDIF
            ENDIF

            IF (i1.EQ.i2) CYCLE

            DO v3 = 1, 3
              v4 = v3 + 1        
              IF (v3.EQ.3) v4 = 1


              IF (i1.EQ.11590.AND.v1.EQ.3) THEN
                dist1 = DSQRT( (ver(tri(i1)%ver(v1),1) - 
     .                          ver(tri(i2)%ver(v3),1))**2 + 
     .                         (ver(tri(i1)%ver(v1),2) - 
     .                          ver(tri(i2)%ver(v3),2))**2)
                dist2 = DSQRT( (ver(tri(i1)%ver(v2),1) - 
     .                          ver(tri(i2)%ver(v4),1))**2 + 
     .                         (ver(tri(i1)%ver(v2),2) - 
     .                          ver(tri(i2)%ver(v4),2))**2)
                dist3 = DSQRT( (ver(tri(i1)%ver(v1),1) - 
     .                          ver(tri(i2)%ver(v4),1))**2 + 
     .                         (ver(tri(i1)%ver(v1),2) - 
     .                          ver(tri(i2)%ver(v4),2))**2)
                dist4 = DSQRT( (ver(tri(i1)%ver(v2),1) - 
     .                          ver(tri(i2)%ver(v3),1))**2 + 
     .                         (ver(tri(i1)%ver(v2),2) - 
     .                          ver(tri(i2)%ver(v3),2))**2)
                IF (dist1.LT.0.0001D0.OR.dist2.LT.0.0001D0.OR.
     .              dist3.LT.0.0001D0.OR.dist4.LT.0.0001D0) THEN
                  WRITE(88,*) 'DIST:',dist1,dist2,i2
                  WRITE(88,*) '    :',dist3,dist4
                  WRITE(88,*) '         ',ver(tri(i1)%ver(v1),1),
     .                                    ver(tri(i1)%ver(v2),1)
                ENDIF
              ENDIF

              IF ((tri(i1)%ver(v1).EQ.tri(i2)%ver(v3).AND.
     .             tri(i1)%ver(v2).EQ.tri(i2)%ver(v4)).OR.
     .            (tri(i1)%ver(v1).EQ.tri(i2)%ver(v4).AND.
     .             tri(i1)%ver(v2).EQ.tri(i2)%ver(v3))) THEN

c...            Adding mapping:
                tri(i1)%map(v1) = i2
                tri(i1)%sid(v1) = v3
                tri(i2)%map(v3) = i1
                tri(i2)%sid(v3) = v1

c...            Neighbour has been found, trigger exit condtion:
                iscan = 999
                EXIT

              ENDIF
            ENDDO

          ENDDO

        ENDDO
      ENDDO


      WRITE(eirfp,*) '  MAPPING SIDES TO SURFACES' 

c...  Map triangles to surfaces:     
      DO i1 = 1, ntri
        IF (tri(i1)%index(1).EQ.1) 
     .    WRITE(eirfp,*) 'RING=',tri(i1)%index(1:2)
        tri(i1)%sur(1:3) = 0
        DO v1 = 1, 3
c          surface_assigned = .FALSE.
          wall_assignment = .FALSE.
          grid_assignment = .FALSE.

          v2 = v1 + 1
          IF (v1.EQ.3) v2 = 1          

          tri(i1)%sur(v1) = 0

          DO i2 = 1, nsurface            
            IF     (surface(i2)%type.EQ.NON_DEFAULT_STANDARD) THEN

              IF (tri(i1)%type.EQ.MAGNETIC_GRID) THEN

                knot   = tri(i1)%index(1)  ! Parameter ...?
                ring   = tri(i1)%index(2)  ! Parameter ...?
                side   = tri(i1)%sideindex(1,v1)
                target = tri(i1)%sideindex(2,v1)

                IF     (surface(i2)%subtype .EQ.STRATUM) THEN
                  IF (surface(i2)%index(1).LE.ring.AND.
     .                surface(i2)%index(2).GE.ring.AND.
     .                surface(i2)%index(3).EQ.target) THEN
                    IF (wall_assignment) THEN
                      CALL WN('ProcessTriangles_06','Wall '//
     .                        'assignment for fluid cell')
                      WRITE(0,*) '  STRATUM:',i1,tri(i1)%index(1:2)
                      STOP  !  Added this on 07/03/2011 -- why should this happen? -SL
                    ELSE
c                      Removed on 07/03/2011 and should be case dependent, i.e. if there are 
c                      triangles behind the target then each particular situation needs to be
c                      processed appropriated. -SL
c                      tri(i1)%map(v1) = 0 ! This should only be set to 0 if the target is opaque...
c                      tri(i1)%sid(v1) = 0 ! ditto
                      tri(i1)%sur(v1) = surface(i2)%num
                      grid_assignment = .TRUE.
                    ENDIF
                  ENDIF
                ELSEIF (surface(i2)%subtype.EQ.
     .                  MAGNETIC_GRID_BOUNDARY) THEN
                    IF (surface(i2)%index(1).LE.knot.AND.
     .                  surface(i2)%index(2).GE.knot.AND.
     .                  surface(i2)%index(3).EQ.ring.AND.
     .                  surface(i2)%index(4).EQ.side) THEN
                      IF (wall_assignment) THEN
                        CALL WN('ProcessTriangles_06','Wall '//
     .                          'assignment for fluid cell')
                        WRITE(0,*) '  BOUNDARY:',i1,tri(i1)%index(1:2)
                      ELSE
                        IF (surface(i2)%index(10).NE.0) THEN
c...         *** MASSIVE HACK : REMAPPING ***
                          DO i3 = 1, nsurface
                            IF (surface(i3)%type.NE.VESSEL_WALL.OR.
     .                          surface(i3)%index(1 ).NE.
     .                          surface(i2)%index(10)) CYCLE
                            tri(i1)%sur(v1) = surface(i3)%index(3)
                            tri(i1)%sideindex(3,v1)=surface(i3)%index(1)
                            tri(i1)%sideindex(4,v1)=surface(i3)%index(2)
                            tri(i1)%sideindex(5,v1)=surface(i3)%index(1)
                            EXIT
                          ENDDO
                          WRITE(0,*) '*** SUPERHACK ***',i3,
     .                      surface(i3)%index(3),
     .                      surface(surface(i3)%index(3))%num
                          IF (i3.EQ.nsurface+1)
     .                      CALL ER('ProcessTriangles_06','Massive re'//
     .                              '-mapping hack failed',*99)
                        ELSE
                          tri(i1)%sur(v1) = surface(i2)%num 
                          grid_assignment = .TRUE.
                          EXIT
                        ENDIF
                      ENDIF
                    ENDIF
                ELSEIF (surface(i2)%subtype .EQ.ADDITIONAL) THEN
c                 Do nothing, these represent the wall in EIRENE:
                ELSE
                  CALL ER('ProcessTriangles_06','Unknown subtype',*99)
                ENDIF
              ENDIF

            ELSEIF (surface(i2)%type.EQ.VESSEL_WALL) THEN
              test = .TRUE.  ! ***REMOVE***
c...          Assign surface end points:
              x(0) = surface(i2)%v(1,1)
              y(0) = surface(i2)%v(2,1)
              x(1) = surface(i2)%v(1,2)
              y(1) = surface(i2)%v(2,2)
c...          Side vertex 1:
              x(2) = ver(tri(i1)%ver(v1),1)  
              y(2) = ver(tri(i1)%ver(v1),2)
              test = test.AND.PointOnLine(x,y,s,t,6,output)
              IF (i1.EQ.problem_triangle) THEN
                WRITE(eirfp,*) '=1>',x(0),y(0)
                WRITE(eirfp,*) '   ',x(1),y(1)
                WRITE(eirfp,*) '   ',x(2),y(2)
                WRITE(eirfp,*) '   ',s,t,test
                IF (test) WRITE(eirfp,*) '     *** MATCH ***'
              ENDIF
c...          Side vertex 2:
              IF (test) THEN
                x(2) = ver(tri(i1)%ver(v2),1)
                y(2) = ver(tri(i1)%ver(v2),2)
                test = test.AND.PointOnLine(x,y,s,t,6,output)
                IF (i1.EQ.problem_triangle) THEN
                  WRITE(eirfp,*) '=2>',x(0),y(0)
                  WRITE(eirfp,*) '   ',x(1),y(1)
                  WRITE(eirfp,*) '   ',x(2),y(2)
                  WRITE(eirfp,*) '   ',s,t,test
                  WRITE(eirfp,*) '   ',i1,v1
                  dummy_test = PointOnLine(x,y,s,t,6,.TRUE.)
                ENDIF
c...            Assign surface index, as appropriate:
                IF (test) THEN
                  IF (i1.EQ.problem_triangle) WRITE(eirfp,*) '   DONE'
                  IF (wall_assignment.OR.grid_assignment) THEN
                    CALL WN('ProcessTriangles_06','More than one '//
     .                      'surface assignment identified C')
                    WRITE(eirfp,*) 'multi-assignment',i1,v1
                    WRITE(eirfp,*) 'current  ',surface(i2)%index(3)
                    WRITE(eirfp,*) '         ',tri(i1)%sideindex(1:3,v1)
                    WRITE(eirfp,*) 'new      ',tri(i1)%sur(v1)
                    WRITE(eirfp,*) '         ',surface(i2)%index(1:2)
                    CYCLE
                  ENDIF
c     .              CALL ER('ProcessTriangles_06','More than one '//
c     .                      'surface assignment identified C',*98)
                  wall_assignment = .TRUE.
                  IF (tri(i1)%map(v1).EQ.-1) THEN
                    tri(i1)%map(v1) = 0
                    tri(i1)%sid(v1) = 0
                  ENDIF
c...              Need to identify which non-default standard surface
c                 should be defined: ... 
                  IF (surface(i2)%index(3).NE.0) THEN
                    tri(i1)%sur(v1) = surface(i2)%index(3)
                    tri(i1)%sideindex(3,v1) = surface(i2)%index(1)  ! Store xVESM index of surface
                    tri(i1)%sideindex(4,v1) = surface(i2)%index(2)  ! Store additional surface index 
                    tri(i1)%sideindex(5,v1) = surface(i2)%index(1)  ! And again, but now the same for all triangles, not just those outside the grid...
                  ELSE
                    CALL ER('ProcessTriangles_06','Surface index '//
     .                      'not assigned',*99)
                  ENDIF
                ENDIF
              ENDIF

            ELSEIF (surface(i2)%type.EQ.HOLE_IN_GRID) THEN
c             Do nothing:
            ELSE
              CALL ER('ProcessTriangles_06','Invalid surface type',*99)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      WRITE(eirfp,*) '  CHECKING ASSIGNMENTS',ntri 



      DO i1 = 1, ntri
        DO v1 = 1, 3                             ! Perhaps eliminate maps through solid surfaces? 
c...      Identify which surface this triangle side is associated with, if any:
          DO i2 = 1, nsurface
            IF (surface(i2)%type    .EQ.NON_DEFAULT_STANDARD.AND.
     .          surface(i2)%index(5).EQ.tri(i1)%sur(v1)) EXIT
          ENDDO       
          IF (i2.EQ.nsurface+1) i2 = 0

          IF (tri(i1)%map(v1).EQ.-1) THEN 
c...        If this shows up again, it may be related to DTOL in PointOnLine:
            WRITE(eirfp,*) 'PROBLEMS WITH MAP (1)',i1  ! NEW EDIT !
            CALL WriteEireneTriangles
            CALL SaveTriangles_06
            CALL DumpGrid('PROBLEM #2 WITH TRIANGLE MAP')

          ELSEIF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.            
     .            tri(i1)%sideindex(2,v1).GT.0.AND.   ! target index
     .            tri(i1)%map      (  v1).NE.0.AND.
     .            tri(i1)%sur      (  v1).NE.0) THEN

c           A triangle side along a target is mapped to a triangle behind the
c           target.  At present, this behind-the-target triangle will not be
c           mapped to the target surface, so make the correction here: 

            i3 = tri(i1)%map(v1)  ! triangle mapped to
            i4 = tri(i1)%sid(v1)  ! and the side of the that triangle

            IF (tri(i3)%sur(i4).EQ.0) THEN
                        
              WRITE(0,*) 'GO:',tri(i1)%map(v1),tri(i1)%sur(v1),
     .                         tri(i1)%index(1:2),
     .                         tri(i1)%type.EQ.MAGNETIC_GRID

              tri(i3)%sur(i4) =  tri(i1)%sur(v1)

            ELSE
              CALL ER('ProcessTriangles_06','Target mapping not '//
     .                'understood',*99)
            ENDIF

          ELSEIF (tri(i1)%map(v1).NE.0.AND.
     .            tri(i1)%sur(v1).NE.0.AND.
     .            i2.GT.0.AND.
     .            surface(MAX(1,i2))%iliin.NE.-1) THEN  ! ILIIN.EQ.-1 = transparent surface
c...        This is being disallowed for now, although there's no
c           reason not to let it happen.  The concern is that, at present, it's
c           more likely to be the result of an error in triangle side to
c           surface mapping than an actual design feature: - SL, 27.07.09
c            CALL ER('ProcessTriangles_06','Map extends through a '//
c     .              'solid surface',*97)
c            WRITE(0,*) 'WARNING ProcessTriangles_06: Map extends '//
c     .                  'through a solid surface, removing link',i3,v3
c            i3 = tri(i1)%map(v1)
c            v3 = tri(i1)%sid(v1)
c            tri(i1)%map(v1) = 0
c            tri(i1)%sid(v1) = 0
c            tri(i3)%map(v3) = 0
c            tri(i3)%sid(v3) = 0

            IF (.NOT.warning_given) THEN
              WRITE(0,*) 'WARNING ProcessTriangles_06: Map extends '//
     .                    'through a solid surface, doing nothing'
              warning_given = .TRUE.
            ENDIF
            WRITE(eirfp,*) 'WARNING ProcessTriangles_06: Map extends'//
     .                     'through solid surface, doing nothing',i1,v1

          ELSEIF (tri(i1)%map(v1).EQ.0) THEN
            IF ( tri(i1)%sur(v1).EQ.0.OR.
     .          (tri(i1)%sur(v1).NE.0.AND.
     .           surface(MAX(1,tri(i1)%sur(v1)))%iliin.EQ.-1.AND.
     .           eir_pass.GT.1)) THEN
c            IF (tri(i1)%map(v1).EQ.-1) THEN 
c...          If this shows up again, it may be related to DTOL in PointOnLine:
c             Indeed, it showed up again, with the ITER grid iterm.carre.105, see above note
c             related to malformed cells. -SL, 23.07.09
              WRITE(eirfp,*) 'PROBLEMS WITH MAP (2)',i1,v1  ! NEW EDIT !
              WRITE(eirfp,*) ' >',tri(i1)%sur(1:3)
              WRITE(eirfp,*) ' >',surface(MAX(1,tri(i1)%sur(1:3)))%iliin
              WRITE(eirfp,*) ' >',tri(i1)%map(1:3)
              WRITE(eirfp,*) ' >',tri(i1)%sid(1:3)

              WRITE(eirfp,*) ' > index     ',tri(i1)%index(1:2)
              WRITE(eirfp,*) ' > sideindex ',tri(i1)%sideindex(1:2,v1)

              WRITE(eirfp,*) '->',ver(tri(i1)%ver(1),1:2)
              WRITE(eirfp,*) ' >',ver(tri(i1)%ver(2),1:2)
              WRITE(eirfp,*) ' >',ver(tri(i1)%ver(3),1:2)
c              WRITE(eirfp,*) '->',ver(tri(380)%ver(1),1:2)
c              WRITE(eirfp,*) ' >',ver(tri(380)%ver(2),1:2)
c              WRITE(eirfp,*) ' >',ver(tri(380)%ver(3),1:2)

              DO i2 = 1, nsurface
                WRITE(eirfp,*) ' ----------'
                WRITE(eirfp,*) ' > subtype ',surface(i2)%subtype,
     .                             STRATUM,MAGNETIC_GRID_BOUNDARY
                WRITE(eirfp,*) ' >   index ',surface(i2)%index(1:4)
                WRITE(eirfp,*) ' >   index ',TRIM(surface(i2)%surtxt)
              ENDDO

              CALL WriteEireneTriangles
              CALL SaveTriangles_06
              CALL DumpGrid('PROBLEM #1 WITH TRIANGLE MAP')
            ENDIF
          ENDIF
c...      Check if 2 close lying points got merged by accident:
          v2 = v2 + 1
          IF (v1.EQ.3) v2 = 1
          IF (tri(i1)%ver(v1).EQ.tri(i1)%ver(v2)) 
     .      CALL ER('ProcessTriangles_06','2 sided triangle',*98)

        ENDDO
      ENDDO

      WRITE(eirfp,*) 'DONE'

      DEALLOCATE(xregion) ! Move to triangle construct since this could be useful elsewhere...?
      DEALLOCATE(yregion)
      DEALLOCATE(nregion)
      DEALLOCATE(iregion)

      RETURN
 97   WRITE(0,'(3X,3I6,2X,I6,2X,A)')
     .  i1,v1,i2,surface(i2)%iliin,TRIM(surface(i2)%surtxt)
      CALL WriteEireneTriangles
      CALL SaveTriangles_06
      CALL DumpGrid('SOLID SURFACE MAP')
      STOP
 98   WRITE(0,*) '  TRIANGLE:',i1,v1
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE OutputVoidSpecification(fp,opt)
      USE mod_sol28
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER, INTENT(IN) :: fp
      TYPE(type_options_eirene), INTENT(IN) :: opt

      INTEGER ivoid

      WRITE(fp,*)
      WRITE(fp,'(A,F10.2)') 'EIRENE TRIANGLE GRID SETTINGS:',
     .  opt%void_version
      IF (opt%void_version.EQ.1.0) THEN 
        WRITE(fp,'(A4,3(2X,A8),A10,2X,A18,2X,A4,2X,3A10)') 
     .    'zone','grid1,2','wall1,2','add1,2','res','hole','code',
     .    'ne','Te','Ti'
        DO ivoid = 1, opt%nvoid
          WRITE(fp,'(I4,3(2X,2I4),F10.4,2X,2F9.4,2X,I4,2X,1P,
     .               3E10.2,0P)') 
     .      opt%void_zone(    ivoid),
     .      opt%void_grid(1:2,ivoid),
     .      opt%void_wall(1:2,ivoid),
     .      opt%void_add (1:2,ivoid),
     .      opt%void_res (    ivoid),
     .      opt%void_hole(1:2,ivoid),
     .      opt%void_code(    ivoid),
     .      opt%void_ne  (    ivoid),
     .      opt%void_te  (    ivoid),
     .      opt%void_ti  (    ivoid)
        ENDDO
      ELSE
        WRITE(fp,'(A4,3(2X,A10),5(2X,A10))') 
     .    'zone','grid','wall','add','res','code',
     .    'ne','Te','Ti'
        DO ivoid = 1, opt%nvoid
          WRITE(fp,'(I4,3(2X,A10),2X,F10.4,2X,I10,1P,3(2X,E10.2),0P)') 
     .      opt%void_zone (ivoid),
     .      TRIM(opt%void2_grid(ivoid)),
     .      TRIM(opt%void2_wall(ivoid)),
     .      TRIM(opt%void2_add (ivoid)),
     .      opt%void_res  (ivoid),
     .      opt%void_code (ivoid),
     .      opt%void_ne   (ivoid),
     .      opt%void_te   (ivoid),
     .      opt%void_ti   (ivoid)
        ENDDO
      ENDIF

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      LOGICAL FUNCTION PointInVoid(x1,y1,nseg,seg,MAXNSEG,
     .                             npts,pts,MAXNPTS,debug,fp)
      IMPLICIT none

      INTEGER, INTENT(IN) :: MAXNPTS,MAXNSEG,nseg,npts,seg(0:MAXNSEG,4),
     .                       fp
      LOGICAL, INTENT(IN) :: debug
      REAL*8 , INTENT(IN) :: x1,y1,pts(MAXNPTS,2)

      REAL*8 , PARAMETER :: DTOL=1.0D-07

      INTEGER iseg,ninter
      REAL*8  x2,x3,x4,y2,y3,y4,s12,s34

      PointInVoid = .FALSE.

      x2 = x1 + 20.0D0
      y2 = y1
      ninter = 0
      DO iseg = 1, nseg
        x3 = pts(seg(iseg,1),1)
        y3 = pts(seg(iseg,1),2)
        x4 = pts(seg(iseg,2),1)
        y4 = pts(seg(iseg,2),2)
        CALL CalcInter(x1,y1,x2,y2,x3,y3,x4,y4,s12,s34) 
        IF (s12.GT.DTOL.AND.s34.GT.0.0D0.AND.s34.LT.1.0D0) 
     .    ninter = ninter + 1
        IF (debug) WRITE(fp,'(4X,A,2F14.7,I4,2F12.5)')
     .      '           :',s12,s34,ninter,x1,y1
      ENDDO  

      IF (ninter.GT.0.AND.MOD(ninter+1,2).EQ.0) PointInVoid = .TRUE.

      RETURN
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE SetupVoidProcessing(opt)
      USE mod_sol28
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      TYPE(type_options_eirene), INTENT(INOUT) :: opt

      INTEGER fp,ivoid,isrf,i1,i2,idefault
      LOGICAL void_specified

      TYPE(type_options_eirene) :: opt_tmp

c...  Output:
      CALL OutputVoidSpecification(eirfp,opt)

      CALL ProcessVoid(-999,opt)  ! Initialisation

c      default_index = 0     
c      DO ivoid = 1, opt%nvoid
c        IF (opt%void_zone.EQ.-1) default_index = ivoid
c      ENDDO

c...  Assign defaults, if necessary:
      idefault  = 0
      void_specified = .FALSE.
      DO i1 = 1, opt%nvoid
        IF (opt%void_zone(i1).EQ.-1) idefault = i1
        IF (opt%void_zone(i1).GE. 1) void_specified = .TRUE.
      ENDDO
      IF (void_specified.AND.idefault.NE.0) 
     .  CALL ER('SetupVoidProcessing','Cannot specify default void '//
     .          'setup and also set specific void parameters',*99)
      IF (idefault.NE.0) THEN
        i1 = idefault
        opt_tmp = opt

        DO isrf = 1, nsurface

          IF (opt%void_version.EQ.1.0) THEN 
            IF (surface(isrf)%type    .NE.NON_DEFAULT_STANDARD  .OR.
     .          surface(isrf)%subtype .NE.MAGNETIC_GRID_BOUNDARY.OR.
     .          surface(isrf)%index(6).LT.opt_tmp%void_grid(1,i1)) CYCLE

            opt%nvoid = opt%nvoid + 1
            i2 = opt%nvoid

            opt%void_grid(:,i2) = surface(isrf)%index(6)
            opt%void_wall(:,i2) = opt_tmp%void_wall(:,i1)
            opt%void_add (:,i2) = opt_tmp%void_add (:,i1)
          ELSE
            IF (surface(isrf)%type    .NE.NON_DEFAULT_STANDARD  .OR.
     .          surface(isrf)%subtype .NE.MAGNETIC_GRID_BOUNDARY.OR.
     .          surface(isrf)%index(6).LE.0) CYCLE

            opt%nvoid = opt%nvoid + 1
            i2 = opt%nvoid

            WRITE(opt%void2_grid(i2),'(I2)') surface(isrf)%index(6)
            opt%void2_wall(i2) = opt_tmp%void2_wall(i1)
            opt%void2_add (i2) = opt_tmp%void2_add (i1)
          ENDIF

          opt%void_zone(i2) = i2 - opt_tmp%nvoid 
          opt%void_res (  i2) = opt_tmp%void_res (  i1)
          opt%void_hole(:,i2) = opt_tmp%void_hole(:,i1)
          opt%void_code(  i2) = opt_tmp%void_code(  i1)
          opt%void_ne  (  i2) = opt_tmp%void_ne  (  i1)
          opt%void_te  (  i2) = opt_tmp%void_te  (  i1)
          opt%void_ti  (  i2) = opt_tmp%void_ti  (  i1)
          WRITE(opt%void_tag(i2),'(A)') 'system default'

        ENDDO

      ENDIF

c...  Output:
      CALL OutputVoidSpecification(eirfp,opt)

c...  Check void blacklist:
      DO i1 = 1, opt%nvoid
        IF (opt%void_zone(i1).NE.-2) CYCLE
        DO i2 = 1, opt%nvoid
          IF (i1.EQ.i2) CYCLE
          IF (opt%void_grid(1,i1).LE.opt%void_zone(i2).AND.
     .        opt%void_grid(2,i1).GE.opt%void_zone(i2)) 
     .      opt%void_zone(i2) = -999  ! Delete
        ENDDO
      ENDDO

c...  Output:
      CALL OutputVoidSpecification(eirfp,opt)

c...  Delete tagged void regions:
      DO i1 = opt%nvoid, 1, -1
        IF (opt%void_zone(i1).NE.-999) CYCLE
        DO i2 = i1, opt%nvoid-1
          opt%void_zone(  i2) = opt%void_zone(  i2+1)-1
          opt%void_grid(:,i2) = opt%void_grid(:,i2+1)
          opt%void_wall(:,i2) = opt%void_wall(:,i2+1)
          opt%void_add (:,i2) = opt%void_add (:,i2+1)
          opt%void_res (  i2) = opt%void_res (  i2+1)
          opt%void_hole(:,i2) = opt%void_hole(:,i2+1)
          opt%void_code(  i2) = opt%void_code(  i2+1)
          opt%void_ne  (  i2) = opt%void_ne  (  i2+1)
          opt%void_te  (  i2) = opt%void_te  (  i2+1)
          opt%void_ti  (  i2) = opt%void_ti  (  i2+1)
        ENDDO
        opt%nvoid = opt%nvoid - 1
      ENDDO

c...  Output:
      CALL OutputVoidSpecification(eirfp,opt)


      RETURN
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE ProcessVoid(izone,opt)
      USE mod_sol28
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER, INTENT(IN) :: izone
      TYPE(type_options_eirene), INTENT(INOUT) :: opt

      LOGICAL PointInVoid,CheckIndex

      INTEGER, PARAMETER :: MAXNSEG = 10000, MAXNPTS = 20000,
     .                      IKLO    = 1    , IKHI    = 2    ,
     .                      MAXNSRF = 10000  ! gfortran, 64-bit
      REAL*8 , PARAMETER :: DTOL=1.0D-06

      INTEGER   fp,ivoid,isrf,isrf1,isrf2,i1,i2,itri,v1,v2,code,
     .          nseg,seg(0:MAXNSEG,7),icnt,nhole,npts,pass,
     .          i3,i4,i5,iseg1,iseg2,ilink,tmp_nseg,ninter,icell,
     .          index,index1,index2,zone_list(100),zone_n
      LOGICAL   debug,cont,link
      CHARACTER command*512,range*128
      REAL      area,ne,te,ti,res(MAXNSRF)
      REAL*8    x1,x2,y1,y2,len,t,tstep,xhole(50),yhole(50),
     .          pts(MAXNPTS,2)

      SAVE

      debug = .TRUE.
      fp = 88



      IF (debug) WRITE(fp,*) 'HERE IN PROCESS VOID',izone

      IF (izone.LT.0) THEN

        IF (nsurface.GT.MAXNSRF)
     .    CALL ER('ProcessVoid_v1_0','Increase MAXNSRF')

        res = 0.0  ! Initialisation

        zone_n    = 0
        zone_list = 0

c...    Loop over the wall surfaces to make sure they match up exactly 
c       with the target segments:
        DO itri = 1, ntri
          DO v1 = 1, 3
            x1 = 0.0D0
            IF (tri(itri)%sideindex(2,v1).EQ.IKLO.OR.
     .          tri(itri)%sideindex(2,v1).EQ.IKHI) THEN
              v2 = v1 + 1
              IF (v1.EQ.3) v2 = 1            
              x1 = ver(tri(itri)%ver(v1),1) 
              y1 = ver(tri(itri)%ver(v1),2)
              x2 = ver(tri(itri)%ver(v2),1)
              y2 = ver(tri(itri)%ver(v2),2)
            ENDIF
            IF (x1.NE.0.0D0) THEN
c...          Search wall surfaces for matching vertices and set them 
c             equal to the grid values, to make sure there are no gaps
c             in the standard wall:
              DO isrf = 1, nsurface
                IF (surface(isrf)%type.NE.VESSEL_WALL) CYCLE
                IF     (DABS(surface(isrf)%v(1,1)-x1).LT.DTOL.AND.
     .                  DABS(surface(isrf)%v(2,1)-y1).LT.DTOL) THEN
                  surface(isrf)%v(1,1) = x1
                  surface(isrf)%v(2,1) = y1
                ELSEIF (DABS(surface(isrf)%v(1,2)-x1).LT.DTOL.AND.
     .                  DABS(surface(isrf)%v(2,2)-y1).LT.DTOL) THEN
                  surface(isrf)%v(1,2) = x1
                  surface(isrf)%v(2,2) = y1
                ENDIF
                IF     (DABS(surface(isrf)%v(1,1)-x2).LT.DTOL.AND.
     .                  DABS(surface(isrf)%v(2,1)-y2).LT.DTOL) THEN
                  surface(isrf)%v(1,1) = x2
                  surface(isrf)%v(2,1) = y2
                ELSEIF (DABS(surface(isrf)%v(1,2)-x2).LT.DTOL.AND.
     .                  DABS(surface(isrf)%v(2,2)-y2).LT.DTOL) THEN
                  surface(isrf)%v(1,2) = x2
                  surface(isrf)%v(2,2) = y2
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDDO 

        RETURN
      ENDIF

c...  Make sure this zone hasn't been process already:
      DO i1 = 1, zone_n
        IF (zone_list(i1).EQ.izone) THEN 
          WRITE(0,*) 'Not processing zone again...',izone
          RETURN
        ENDIF
      ENDDO
      zone_n = zone_n + 1
      zone_list(zone_n) = izone

      nseg  = 0
      seg   = 0
      npts  = 0       
      nhole = 0
      area  = 0.0
      ne    = 0.0
      te    = 0.0
      ti    = 0.0

      DO ivoid = 1, opt%nvoid
        IF (opt%void_zone(ivoid).NE.izone) CYCLE

        IF (debug) WRITE(fp,*) '  PROCESSING VOID SETUP',ivoid
c       ----------------------------------------------------------------
c...    Examine the outer radial boundary surfaces of the fluid grid and 
c       collect the associated line segments:

c       The boundary surface indices have been identified, so search the
c       list of triangles for sides that match up with these surfaces:

        range = TRIM(opt%void2_grid(ivoid))

        IF (range(1:4).NE.'none') THEN

          DO itri = 1, ntri
            DO v1 = 1, 3
              isrf = tri(itri)%sur(v1)
              IF (isrf.EQ.0) CYCLE

              index = surface(isrf)%index(6)  

              WRITE(0,*) 'range ',range
              WRITE(0,*) 'index ',index

              IF (surface(isrf)%type    .EQ.NON_DEFAULT_STANDARD  .AND.
     .            (surface(isrf)%subtype.EQ.MAGNETIC_GRID_BOUNDARY.OR.
     .            (surface(isrf)%subtype.EQ.STRATUM)).AND.              
     .            CheckIndex(index,range)) THEN

                nseg = nseg + 1
                seg(nseg,1) = npts + 1
                seg(nseg,2) = npts + 2
                seg(nseg,5) = index
                v2 = v1 + 1
                IF (v1.EQ.3) v2 = 1
                npts = npts + 1
                pts(npts,1) = ver(tri(itri)%ver(v1),1) 
                pts(npts,2) = ver(tri(itri)%ver(v1),2)
                npts = npts + 1
                pts(npts,1) = ver(tri(itri)%ver(v2),1)
                pts(npts,2) = ver(tri(itri)%ver(v2),2)
                IF (debug) THEN
                  WRITE(fp,'(A,2I6,2X,A)') 'GRID:',isrf,index,range
                  WRITE(fp,'(5I6)') nseg,itri,isrf,isrf1,isrf2
                  WRITE(fp,*) '    PTS1=',pts(npts-1,1:2)
                  WRITE(fp,*) '    PTS2=',pts(npts  ,1:2)
                ENDIF
              ENDIF
            ENDDO              
          ENDDO      

          IF (nseg.EQ.0)  
     .      CALL ER('ProcessVoid','No grid related surfaces found',*98)
        ENDIF


c       ------------------------------------------------------------------
c...    Search through the list of standard wall line segments and build
c       a list for each zone that completes the individual voids:
        range = TRIM(opt%void2_wall(ivoid))

        WRITE(fp,*) 'WALL I1,2-= ',i1,i2

        IF (range(1:3).EQ.'def') THEN
c         Find the wall segments for this zone automatically.  Just keep
c         mindlessly filing through the wall segments until the path is 
c         closed:
          cont = .TRUE.
          pass = 0
          DO WHILE(cont)
            pass = pass + 1
            IF (pass.EQ.100) STOP 'NOT PASSING...'
            cont = .FALSE.
c           Check if there is already a link to this segment:
            tmp_nseg = nseg
            DO iseg1 = 1, tmp_nseg
              IF (debug) WRITE(fp,*) '  Trying1-:',iseg1,seg(iseg1,3)
              IF (seg(iseg1,3).EQ.1) CYCLE
c             Check both ends of the current focus segment:
              DO ilink = 1, 2
                link = .FALSE.
                DO iseg2 = 1, nseg
                  IF (debug) WRITE(fp,*) '  Trying2-:',iseg2
                  IF (iseg1.EQ.iseg2) CYCLE
                  i3 = seg(iseg1,ilink)
c                 Check both ends of the test segment:
                  i4 = seg(iseg2,1)
                  i5 = seg(iseg2,2)
                  IF ((DABS(pts(i3,1)-pts(i4,1)).LT.DTOL.AND.
     .                 DABS(pts(i3,2)-pts(i4,2)).LT.DTOL).OR.
     .                (DABS(pts(i3,1)-pts(i5,1)).LT.DTOL.AND.
     .                 DABS(pts(i3,2)-pts(i5,2)).LT.DTOL)) THEN
                    link = .TRUE.
                    IF (ilink.EQ.2) seg(iseg1,3) = 1
                    EXIT
                  ENDIF
                ENDDO
                IF (debug) WRITE(fp,*) '  PASS:',pass,iseg1,link
                IF (.NOT.link) EXIT
              ENDDO

c             Both ends of the segment are attached so start
c             looking at the next segment:
              IF (link) CYCLE

              cont = .TRUE.
              DO isrf = 1, nsurface
                WRITE(fp,*) 'Trying 3-:',isrf,
     .                       surface(isrf)%type,VESSEL_WALL,
     .                       seg(iseg1,4)
                IF (surface(isrf)%type.NE.VESSEL_WALL.OR.    
     .              seg(iseg1,4).EQ.isrf) CYCLE
                x1 = surface(isrf)%v(1,1)
                y1 = surface(isrf)%v(2,1)
                x2 = surface(isrf)%v(1,2)
                y2 = surface(isrf)%v(2,2)
                IF ((DABS(pts(i3,1)-x1).LT.DTOL.AND.
     .               DABS(pts(i3,2)-y1).LT.DTOL).OR.
     .              (DABS(pts(i3,1)-x2).LT.DTOL.AND.
     .               DABS(pts(i3,2)-y2).LT.DTOL)) THEN

                  WRITE(fp,*) 'Wall surface 1-',isrf,
     .                   surface(isrf)%index(1:2)

                  IF (res(isrf).EQ.0.0) THEN
                    res(isrf) = opt%void_res(ivoid)
                  ELSE
                    WRITE(fp,*) 'ISRF,RES:',isrf,nsurface,res(isrf)
                    WRITE(fp,*) 'X1,Y1:',x1,y1
                    WRITE(fp,*) 'X2,Y2:',x2,y2
                    WRITE(fp,*) 'RES  :',res                  
                    STOP 'NOT READY FOR THE RES...'
                  ENDIF

                  len = DSQRT((x1 - x2)**2 + (y1 - y2)**2)
                  IF (len.GT.DBLE(res(isrf))) THEN
                    tstep = 1.0D0 / DBLE(INT(len/DBLE(res(isrf))) + 1)
                  ELSE
                    tstep = 1.0D0
                  ENDIF

                  IF (debug) THEN
                    WRITE(fp,'(A,2I6,2X,2I6,2F10.4)') 
     .                ' NEW WALL SEG:',iseg1,ilink, 
     .                surface(isrf)%index(1:2),res(isrf),tstep

                    WRITE(fp,*) '    I3    =',i3
                    WRITE(fp,*) '    PTS1,2=',pts(i3,1:2)
                    WRITE(fp,*) '    ISEG11=',pts(seg(iseg1,1),1:2)
                    WRITE(fp,*) '    ISEG12=',pts(seg(iseg1,2),1:2)
                    WRITE(fp,*) '    ISRF,RES:',isrf,res(isrf)
                    WRITE(fp,*) '    X1,Y1:',x1,y1
                    WRITE(fp,*) '    X2,Y2:',x2,y2

c                    STOP 'dfsdfsd'
                  ENDIF
              
                  DO t = 0.0D0, 0.9999999D0, tstep 
                    nseg = nseg + 1
                    seg(nseg,1) = npts + 1
                    seg(nseg,2) = npts + 2
                    seg(nseg,4) = isrf
                    seg(nseg,6) = surface(isrf)%index(1)
                    npts = npts + 1
                    pts(npts,1) = x1 + t * (x2 - x1)             ! Orientation is correct
                    pts(npts,2) = y1 + t * (y2 - y1) 
                    npts = npts + 1
                    pts(npts,1) = x1 + (t + tstep) * (x2 - x1) 
                    pts(npts,2) = y1 + (t + tstep) * (y2 - y1)
                  ENDDO
                  EXIT
                ENDIF
              ENDDO
              IF (isrf.EQ.nsurface+1) 
     .          CALL ER('ProcessVoid','No link to wall surface',*99)
            ENDDO

          ENDDO

        ELSEIF (range(1:4).NE.'none') THEN
          DO isrf = 1, nsurface

            index = surface(isrf)%index(1)  

            WRITE(fp,*) 'Trying 4-:',isrf,index,
     .                   surface(isrf)%type,VESSEL_WALL,
     .                   .NOT.CheckIndex(index,range)

            IF (surface(isrf)%type.NE.VESSEL_WALL.OR.
     .          .NOT.CheckIndex(index,range)) CYCLE

            WRITE(fp,*) 'Wall surface 2-',isrf,surface(isrf)%index(1)

            x1 = surface(isrf)%v(1,1)
            y1 = surface(isrf)%v(2,1)
            x2 = surface(isrf)%v(1,2)
            y2 = surface(isrf)%v(2,2)

            IF     (res(isrf).EQ.0.0                ) THEN
              res(isrf) = opt%void_res(ivoid)
            ELSEIF (res(isrf).NE.opt%void_res(ivoid)) THEN
              WRITE(fp,*) 'ISRF,RES:',isrf,nsurface,res(isrf)
              WRITE(fp,*) 'X1,Y1:',x1,y1
              WRITE(fp,*) 'X2,Y2:',x2,y2
              WRITE(fp,*) 'RES  :',res                  
              STOP 'NOT READY FOR THE RES... (2)'
            ENDIF

            len = DSQRT((x1 - x2)**2 + (y1 - y2)**2)
            IF (len.GT.DBLE(res(isrf))) THEN
              tstep = 1.0D0 / DBLE(INT(len/DBLE(res(isrf))) + 1)
            ELSE
              tstep = 1.0D0
            ENDIF

            DO t = 0.0D0, 0.9999999D0, tstep 
              nseg = nseg + 1
              seg(nseg,1) = npts + 1
              seg(nseg,2) = npts + 2
              seg(nseg,4) = isrf
              seg(nseg,6) = surface(isrf)%index(1)
              npts = npts + 1
              pts(npts,1) = x1 + t * (x2 - x1)             ! Orientation is correct
              pts(npts,2) = y1 + t * (y2 - y1) 
              npts = npts + 1
              pts(npts,1) = x1 + (t + tstep) * (x2 - x1) 
              pts(npts,2) = y1 + (t + tstep) * (y2 - y1)
            ENDDO

          ENDDO
        ENDIF
      ENDDO  ! IVOID loop

        DO isrf = 1, nsurface
           WRITE(fp,*) 'res:',isrf,res(isrf),surface(isrf)%index(1:2)
        ENDDO



c     ------------------------------------------------------------------
c...  Check for any additional line segments to include (where 'additional'
c     is from the EIRENE nomenclature meaning extra surfaces not directly
c     related to the fluid grid): 

      DO ivoid = 1, opt%nvoid
        IF (opt%void_zone(ivoid).NE.izone) CYCLE

c        i1 = opt%void_add(1,ivoid)
c        i2 = opt%void_add(2,ivoid)
c        IF (i1.LE.0.OR.i2.LE.0) CYCLE
        range = TRIM(opt%void2_add(ivoid))
 
        IF (range(1:3).EQ.'def'.OR.range(1:4).EQ.'none') CYCLE 

        tmp_nseg = nseg
       
        DO isrf = 1, nsurface
          WRITE(fp,*) 'Add check',surface(isrf)%type.NE.VESSEL_WALL,
     .                surface(isrf)%index(2),i1,i2
        
          index = surface(isrf)%index(2)  

          IF ((surface(isrf)%type.NE.VESSEL_WALL.AND.
     .         surface(isrf)%type.NE.HOLE_IN_GRID).OR.
     .        .NOT.CheckIndex(index,range)) CYCLE
c     .        (surface(isrf)%index(2).LT.i1.OR.
c     .         surface(isrf)%index(2).GT.i2)) CYCLE

          x1 = surface(isrf)%v(1,1)
          y1 = surface(isrf)%v(2,1)
          x2 = surface(isrf)%v(1,2)
          y2 = surface(isrf)%v(2,2)

          IF (surface(isrf)%type.EQ.HOLE_IN_GRID) THEN
            WRITE(fp,*) 'Add hole',isrf,surface(isrf)%index(2)
            nhole = nhole + 1
            xhole(nhole) = x1
            yhole(nhole) = y1
          ELSE
            WRITE(fp,*) 'Add surface',isrf,surface(isrf)%index(2)

c        DO i1 = 1, nsurface
c           WRITE(fp,*) 'res:',i1,res(i1),surface(i1)%index(1:2)
c        ENDDO
            
            IF     (res(isrf).EQ.0.0                ) THEN
              res(isrf) = opt%void_res(ivoid)
            ELSEIF (res(isrf).NE.opt%void_res(ivoid)) THEN
              WRITE(fp,*) 'ISRF,RES:',isrf,nsurface,res(isrf)
              WRITE(fp,*) 'X1,Y1:',x1,y1
              WRITE(fp,*) 'X2,Y2:',x2,y2
              WRITE(fp,*) 'RES  :',res                  
              STOP 'NOT READY FOR THE RES... (3)'
            ENDIF
            
            len = DSQRT((x1 - x2)**2 + (y1 - y2)**2)
            IF (len.GT.DBLE(res(isrf))) THEN
              tstep = 1.0D0 / DBLE(INT(len/DBLE(res(isrf))) + 1)
            ELSE
              tstep = 1.0D0
            ENDIF
            DO t = 0.0D0, 0.9999999D0, tstep 
              nseg = nseg + 1
              seg(nseg,1) = npts + 1
              seg(nseg,2) = npts + 2
              seg(nseg,4) = isrf
              seg(nseg,7) = surface(isrf)%index(2)
              npts = npts + 1
              pts(npts,1) = x1 + t * (x2 - x1)             ! Orientation is correct
              pts(npts,2) = y1 + t * (y2 - y1) 
              npts = npts + 1
              pts(npts,1) = x1 + (t + tstep) * (x2 - x1) 
              pts(npts,2) = y1 + (t + tstep) * (y2 - y1)
            ENDDO
          ENDIF        

        ENDDO  ! ISRF loop

        IF (nseg.EQ.tmp_nseg) 
     .    CALL ER('ProcessVoid','Additional surfaces requested but '//
     .            'none assigned',*98)

      ENDDO  ! IVOID loop


c     ----------------------------------------------------------------------
c...  Sort segments so that the first ??? form a coherent boundary:
      DO i2 = 1, nseg-1
        DO i3 = i2+1, nseg
          IF     (DABS(pts(seg(i2,2),1)-pts(seg(i3,1),1)).LT.DTOL.AND.
     .            DABS(pts(seg(i2,2),2)-pts(seg(i3,1),2)).LT.DTOL) THEN
            IF (i3.EQ.i2+1) THEN
c...          Do nothing, all okay:
              EXIT
            ELSE
              seg(0   ,1:7) = seg(i2+1,1:7)
              seg(i2+1,1:7) = seg(i3  ,1:7)
              seg(i3  ,1:7) = seg(0   ,1:7)
              EXIT
            ENDIF
          ELSEIF (DABS(pts(seg(i2,2),1)-pts(seg(i3,2),1)).LT.DTOL.AND.
     .            DABS(pts(seg(i2,2),2)-pts(seg(i3,2),2)).LT.DTOL) THEN
            IF (i3.EQ.i2+1) THEN
              seg(0 ,1) = seg(i3,1)
              seg(i3,1) = seg(i3,2) ! Swap the order of the points
              seg(i3,2) = seg(0 ,1)
              EXIT
            ELSE
              seg(0   ,1:7) = seg(i2+1,1:7)
              seg(i2+1,1  ) = seg(i3  ,2  )  ! Swap the order of the points
              seg(i2+1,2  ) = seg(i3  ,1  )
              seg(i2+1,3:7) = seg(i3  ,3:7)
              seg(i3  ,1:7) = seg(0   ,1:7)
              EXIT
            ENDIF
          ENDIF
        ENDDO
        IF (i3.EQ.nseg+1) THEN
          WRITE(fp,*) 'trouble at',i2 
          DO i4 = 1, nseg
            WRITE(fp,'(A,I6,2(2X,2F12.6),2X,3I6)') 
     .        'segments ',i4,pts(seg(i4,1),1:2),pts(seg(i4,2),1:2),
     .        seg(i4,5:7)
          ENDDO        
          CALL ER('ProcessVoid','Zone perimeter gap detected',*99)
        ENDIF

c       Check if the path has been closed, assuming that the first segment is a
c       valid path segment (need to keep this in mind!):
c        WRITE(fp,*) 'i2,3:',i2,i3
c        WRITE(fp,*) pts(seg(i2,1),1:2)
c        WRITE(fp,*) pts(seg(i2,2),1:2)
c        WRITE(fp,*) pts(seg(i2+1,1),1:2)
c        WRITE(fp,*) pts(seg(i2+1,2),1:2)
c        DO i4 = 1, nseg
c          WRITE(fp,*) 'segment order',i4,seg(i4,5)   
c        ENDDO

        IF (DABS(pts(seg(i2+1,2),1)-pts(seg(1,1),1)).LT.DTOL.AND.
     .      DABS(pts(seg(i2+1,2),2)-pts(seg(1,1),2)).LT.DTOL) THEN 
          WRITE(fp,*) 'CLOSURE DETECTED!',i2,i3,nseg
          EXIT
        ENDIF
      ENDDO

c...  Register whether or not all of the segments that were loaded were needed
c     to close the path, or whether there are some extra ones:    
      IF (i2.NE.nseg-1) THEN
        WRITE(fp,'(A,I6)') 'termination',i2
        DO i3 = 1, nseg
          WRITE(fp,'(A,I6,2(2X,2F12.6),2X,3I6)') 
     .      'segments ',i3,pts(seg(i3,1),1:2),pts(seg(i3,2),1:2),
     .      seg(i3,5:7)
c          WRITE(fp,*) 'segment order',i2,seg(i3,5)   
        ENDDO
        CALL WN('ProcessVoid','Non standard wall job')
      ENDIF

c      DO i2 = 1, nseg
c        WRITE(fp,*) 'segment order',i2,seg(i2,5)   
c      ENDDO

c     IF (debug) THEN
c        DO i1 = 1, nseg
c          WRITE(fp,*) 'PERIMETER:',i1,pts(seg(i1,1),1:2)
c          WRITE(fp,*) '         :',i1,pts(seg(i1,2),1:2)
c        ENDDO
c      ENDIF

c     ------------------------------------------------------------------
c...  The perimeter of the void / zone is now defined, so look to see if 
c     any additional line segments or holes should be included:     
      tmp_nseg = i2 + 1 ! nseg
      DO ivoid = 1, opt%nvoid
        IF (opt%void_zone(ivoid).NE.izone) CYCLE

c...    Set requested triangle area:
        IF (area.EQ.0.0.AND.opt%void_res(ivoid).NE.0.0) 
     .    area = 0.5 * opt%void_res(ivoid)**2

c...    Plasma conditions for all triangles in zone -- but note that
c       this does not imply any associated recycling in EIRENE:
        IF (opt%void_ne(ivoid).GT.0.0.AND.ne.EQ.0.0) THEN
          ne = opt%void_ne(ivoid)
          te = opt%void_te(ivoid)
          ti = opt%void_ti(ivoid)
        ENDIF


        range = TRIM(opt%void2_add(ivoid))
        IF (range(1:3).EQ.'def') THEN
c        IF (opt%void_add(1,ivoid).EQ.-1.AND.
c     .      opt%void_add(2,ivoid).EQ.-1) THEN

c...      Look for additional surfaces that are inside each zone:
          DO isrf = 1, nsurface
            index = surface(isrf)%index(2)

            IF (surface(isrf)%type    .NE.VESSEL_WALL.OR.
     .          surface(isrf)%index(2).EQ.0          .OR.
     .          res(isrf).NE.0.0) CYCLE

            x1 = surface(isrf)%v(1,1)  ! Only check the first point of each surface,
            y1 = surface(isrf)%v(2,1)  ! assuming the wall isn't ill posed...
            IF (debug) WRITE(fp,'(4X,A,3I4)')
     .        'ADD TRYING:',isrf,surface(isrf)%index(1:2)
            IF (PointInVoid(x1,y1,tmp_nseg,seg,MAXNSEG,
     .                      npts,pts,MAXNPTS,debug,fp)) THEN
c...          The point is inside the void perimeter, so add the segment to the list:
              res(isrf) = opt%void_res(ivoid)
c...          Determine if the segment needs to be sub-divided:
              x1 = surface(isrf)%v(1,1)
              y1 = surface(isrf)%v(2,1)
              x2 = surface(isrf)%v(1,2)
              y2 = surface(isrf)%v(2,2)
              len = DSQRT((x1 - x2)**2 + (y1 - y2)**2)
              IF (len.GT.DBLE(res(isrf))) THEN
                tstep = 1.0D0 / DBLE(INT(len/DBLE(res(isrf))) + 1)
              ELSE
                tstep = 1.0D0
              ENDIF
              IF (debug) THEN
                WRITE(fp,*) ' NEW ADD SEG:',isrf
c                WRITE(fp,*) '    I3    =',i3
c                WRITE(fp,*) '    PTS1,2=',pts(i3,1:2)
c                WRITE(fp,*) '    ISEG11=',pts(seg(iseg1,1),1:2)
c                WRITE(fp,*) '    ISEG12=',pts(seg(iseg1,2),1:2)
c                WRITE(fp,*) '    ISRF,RES:',isrf,res(isrf)
c                WRITE(fp,*) '    X1,Y1:',x1,y1
c                WRITE(fp,*) '    X2,Y2:',x2,y2
c                STOP 'dfsdfsd'
              ENDIF
              DO t = 0.0D0, 0.9999999D0, tstep 
                nseg = nseg + 1
                seg(nseg,1) = npts + 1
                seg(nseg,2) = npts + 2
                seg(nseg,4) = isrf
                npts = npts + 1
                pts(npts,1) = x1 + t * (x2 - x1)             ! Orientation is correct
                pts(npts,2) = y1 + t * (y2 - y1) 
                npts = npts + 1
                pts(npts,1) = x1 + (t + tstep) * (x2 - x1) 
                pts(npts,2) = y1 + (t + tstep) * (y2 - y1)
              ENDDO
            ENDIF
          ENDDO

c...      Look for holes:
          DO isrf = 1, nsurface
            IF (surface(isrf)%type.NE.HOLE_IN_GRID) CYCLE
            x1 = surface(isrf)%v(1,1)  
            y1 = surface(isrf)%v(2,1)  
            IF (debug) WRITE(fp,'(4X,A,3I4)')
     .        'HOLE TRYING:',isrf,surface(isrf)%index(1:2)
            IF (PointInVoid(x1,y1,tmp_nseg,seg,MAXNSEG,
     .                      npts,pts,MAXNPTS,debug,fp)) THEN
              nhole = nhole + 1
              xhole(nhole) = x1
              yhole(nhole) = y1
              IF (debug) WRITE(fp,*) '   ADDING HOLE:',nhole,x1,y1
            ENDIF
          ENDDO

        ENDIF    
 
      ENDDO


c...  Eliminate duplicate verticies:
      DO i2 = 1, npts
        DO i3 = i2+1, npts
          IF (pts(i2,1).NE.-999.0.AND.
     .        DABS(pts(i2,1)-pts(i3,1)).LT.DTOL.AND.
     .        DABS(pts(i2,2)-pts(i3,2)).LT.DTOL) THEN
            pts(i3,1) = -999.0
            pts(i3,2) = -999.0
            DO i4 = 1, nseg
              IF (seg(i4,1).EQ.i3) seg(i4,1) = i2
              IF (seg(i4,2).EQ.i3) seg(i4,2) = i2
            ENDDO
          ENDIF
        ENDDO
      ENDDO
c...  Delete points:
      DO i2 = npts, 1, -1
        IF (pts(i2,1).EQ.-999.0) THEN
          DO i3 = i2, npts-1
            pts(i3,1) = pts(i3+1,1)
            pts(i3,2) = pts(i3+1,2)
          ENDDO
          DO i4 = 1, nseg
            IF (seg(i4,1).GE.i2) seg(i4,1) = seg(i4,1) - 1
            IF (seg(i4,2).GE.i2) seg(i4,2) = seg(i4,2) - 1
          ENDDO
          npts = npts - 1
        ENDIF
      ENDDO


      DO i1 = 1, nseg
        WRITE(fp,'(A,I6,2(2F12.6,2X))') 
     .    'VOID:',i1,pts(seg(i1,1),1),pts(seg(i1,1),2),
     .               pts(seg(i1,2),1),pts(seg(i1,2),2)
      ENDDO


      fp = 99      
      OPEN(UNIT=fp,FILE='triangle.poly',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=99)      
      WRITE(fp,*) npts,2,1,0
      DO i2 = 1, npts
c        WRITE(fp,'(I6,2F12.7)') i2,pts(i2,1),pts(i2,2)
        WRITE(fp,'(I6,2F19.14)') i2,pts(i2,1),pts(i2,2)
      ENDDO
      WRITE(fp,*) nseg,0
      DO i2 = 1, nseg
        IF (i2.EQ.nseg) THEN
          WRITE(fp,'(4I6)') i2,seg(i2,1),seg(i2,2),0
        ELSE
          WRITE(fp,'(4I6)') i2,seg(i2,1),seg(i2,2),1
        ENDIF
      ENDDO
      WRITE(fp,*) nhole
      DO i2 = 1, nhole
        WRITE(fp,'(I6,2F12.7)') i2,xhole(i2),yhole(i2)
      ENDDO
      CLOSE (fp)

c...  Call triangle:
      IF (area.EQ.0.0) area = 0.01
      WRITE(command,10) 'triangle -p -q -a',area,' -Y triangle.poly>tmp'
 10   FORMAT(A,F10.8,A)
      WRITE(eirfp,*) 'COMMAND: >'//command(1:LEN_TRIM(command))//'<'
      WRITE(0    ,*) 'COMMAND: >'//command(1:LEN_TRIM(command))//'<'
      CALL CIssue(command(1:LEN_TRIM(command)),code)
      WRITE(eirfp,*) 'RETURN_CODE:',code

      icnt = icnt + 1

      CALL ReadPolyFile_06(izone,ne,te,ti)


      RETURN
 99   WRITE(0,*) 'I3    =',i3
      WRITE(0,*) 'PTS1,2=',pts(i3,1:2)
 98   CONTINUE
      CALL WriteEireneTriangles
      CALL DumpGrid('PROBLEM WITH VOID')
      STOP
      END
c
c ====================================================================
c
      SUBROUTINE ProcessVoid_v1_0(izone,opt)
      USE mod_sol28
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER, INTENT(IN) :: izone
      TYPE(type_options_eirene), INTENT(INOUT) :: opt

      LOGICAL PointInVoid

      INTEGER, PARAMETER :: MAXNSEG = 10000, MAXNPTS = 20000,
     .                      IKLO    = 1    , IKHI    = 2    ,
     .                      MAXNSRF = 10000  ! gfortran, 64-bit
      REAL*8 , PARAMETER :: DTOL=1.0D-06

      INTEGER   fp,ivoid,isrf,isrf1,isrf2,i1,i2,itri,v1,v2,code,
     .          nseg,seg(0:MAXNSEG,5),icnt,nhole,npts,pass,
     .          i3,i4,i5,iseg1,iseg2,ilink,tmp_nseg,ninter,icell
      LOGICAL   debug,cont,link
      CHARACTER command*512
      REAL      area,ne,te,ti,res(MAXNSRF) !(nsurface)
      REAL*8    x1,x2,y1,y2,len,t,tstep,xhole(50),yhole(50),
     .          pts(MAXNPTS,2)

      SAVE

      debug = .TRUE.
      fp = 88

      IF (debug) WRITE(fp,*) 'HERE IN PROCESS VOID',izone

      IF (izone.LT.0) THEN
        IF (debug) WRITE(fp,*) 'INITIALIZATION'

        IF (nsurface.GT.MAXNSRF)
     .    CALL ER('ProcessVoid_v1_0','Increase MAXNSRF')

        res = 0.0  ! Initialisation

c        DO isrf = 1, nsurface
c           WRITE(fp,*) 'res chk:',isrf,res(isrf),
c     .                            surface(isrf)%index(1:2)
c        ENDDO


c        IF (izone.EQ.-999) res = 0.0  ! Initialisation

c...    Loop over the wall surfaces to make sure they match up exactly 
c       with the target segments:
        DO itri = 1, ntri
          DO v1 = 1, 3
            x1 = 0.0D0
            IF (tri(itri)%sideindex(2,v1).EQ.IKLO.OR.
     .          tri(itri)%sideindex(2,v1).EQ.IKHI) THEN
              v2 = v1 + 1
              IF (v1.EQ.3) v2 = 1            
              x1 = ver(tri(itri)%ver(v1),1) 
              y1 = ver(tri(itri)%ver(v1),2)
              x2 = ver(tri(itri)%ver(v2),1)
              y2 = ver(tri(itri)%ver(v2),2)
            ENDIF
            IF (x1.NE.0.0D0) THEN
c...          Search wall surfaces for matching vertices and set them 
c             equal to the grid values, to make sure there are no gaps
c             in the standard wall:
              DO isrf = 1, nsurface
                IF (surface(isrf)%type.NE.VESSEL_WALL) CYCLE
                IF     (DABS(surface(isrf)%v(1,1)-x1).LT.DTOL.AND.
     .                  DABS(surface(isrf)%v(2,1)-y1).LT.DTOL) THEN
c                  WRITE(0,*) 'WALL FIX 1'
                  surface(isrf)%v(1,1) = x1
                  surface(isrf)%v(2,1) = y1
                ELSEIF (DABS(surface(isrf)%v(1,2)-x1).LT.DTOL.AND.
     .                  DABS(surface(isrf)%v(2,2)-y1).LT.DTOL) THEN
c                  WRITE(0,*) 'WALL FIX 2'
                  surface(isrf)%v(1,2) = x1
                  surface(isrf)%v(2,2) = y1
                ENDIF
                IF     (DABS(surface(isrf)%v(1,1)-x2).LT.DTOL.AND.
     .                  DABS(surface(isrf)%v(2,1)-y2).LT.DTOL) THEN
c                  WRITE(0,*) 'WALL FIX 3'
                  surface(isrf)%v(1,1) = x2
                  surface(isrf)%v(2,1) = y2
                ELSEIF (DABS(surface(isrf)%v(1,2)-x2).LT.DTOL.AND.
     .                  DABS(surface(isrf)%v(2,2)-y2).LT.DTOL) THEN
c                  WRITE(0,*) 'WALL FIX 4'
                  surface(isrf)%v(1,2) = x2
                  surface(isrf)%v(2,2) = y2
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDDO 

        RETURN
      ENDIF

      nseg = 0
      seg  = 0
      npts = 0       
      nhole = 0
      area = 0.0
      ne = 0.0
      te = 0.0
      ti = 0.0

      IF (debug) WRITE(fp,*) 'BEYOND INITIALISATION',izone

      DO isrf = 1, nsurface
         WRITE(fp,*) 'res chk:',isrf,res(isrf),surface(isrf)%index(1:2)
      ENDDO

      DO ivoid = 1, opt%nvoid
        IF (opt%void_zone(ivoid).NE.izone) CYCLE

        IF (debug) WRITE(fp,*) '  PROCESSING VOID SETUP',ivoid
c       ----------------------------------------------------------------
c...    Examine the outer radial boundary surfaces of the fluid grid and 
c       collect the associated line segments:
c
c       Map the surface indices provided in the input file to the surface
c       indices defined in the EIRENE setup code:
        isrf1 = -1
        isrf2 = -1
        i1 = opt%void_grid(1,ivoid)
        i2 = opt%void_grid(2,ivoid)
        DO isrf = 1, nsurface
c          IF (debug) WRITE(fp,*) 'GRID SURFACE:',isrf
          IF (surface(isrf)%type   .NE.NON_DEFAULT_STANDARD  .OR.
     .        surface(isrf)%subtype.NE.MAGNETIC_GRID_BOUNDARY) CYCLE
          IF (debug) WRITE(fp,*) 'GRID SURFACE: OK',isrf,
     .                            surface(isrf)%index(6)
          IF (isrf1.EQ.-1.AND.surface(isrf)%index(6).GE.i1) isrf1 = isrf
          IF (                surface(isrf)%index(6).LE.i2) isrf2 = isrf
        ENDDO

        IF (debug) WRITE(fp,*) 'GRID:',ivoid,i1,i2,isrf1,isrf2

c       The boundary surface indices have been identified, so search the
c       list of triangles for sides that match up with these surfaces:
        IF (isrf1.NE.-1.AND.isrf2.NE.-1) THEN
          DO itri = 1, ntri
            DO v1 = 1, 3
              isrf = tri(itri)%sur(v1)
              IF (isrf.EQ.0) CYCLE
              IF (surface(isrf)%type    .EQ.NON_DEFAULT_STANDARD  .AND.
     .            surface(isrf)%subtype .EQ.MAGNETIC_GRID_BOUNDARY.AND.
     .            surface(isrf)%index(6).GT.0.AND. 
     .            isrf.GE.isrf1.AND.isrf.LE.isrf2) THEN
                nseg = nseg + 1
                seg(nseg,1) = npts + 1
                seg(nseg,2) = npts + 2
                seg(nseg,5) = -1
                v2 = v1 + 1
                IF (v1.EQ.3) v2 = 1
                npts = npts + 1
                pts(npts,1) = ver(tri(itri)%ver(v1),1) 
                pts(npts,2) = ver(tri(itri)%ver(v1),2)
                npts = npts + 1
                pts(npts,1) = ver(tri(itri)%ver(v2),1)
                pts(npts,2) = ver(tri(itri)%ver(v2),2)
                IF (debug) THEN
                  WRITE(fp,'(A,5I6)') 'GRID:',nseg,itri,isrf,isrf1,isrf2
                  WRITE(fp,*) '    PTS1=',pts(npts-1,1:2)
                  WRITE(fp,*) '    PTS2=',pts(npts  ,1:2)
                ENDIF
              ENDIF
            ENDDO              
          ENDDO      
        ENDIF

c       ------------------------------------------------------------------
c...    Search through the list of standard wall line segments and build
c       a list for each zone that completes the individual voids:
        i1 = opt%void_wall(1,ivoid)
        i2 = opt%void_wall(2,ivoid)

        WRITE(fp,*) 'WALL I1,2= ',i1,i2

        IF     (i1.EQ.-1.AND.i2.EQ.-1) THEN
c         Find the wall segments for this zone automatically -- just keep
c         mindlessly filing through the wall segments until the path is 
c         closed:
          cont = .TRUE.
          pass = 0
          DO WHILE(cont)
            pass = pass + 1
            IF (pass.EQ.100) STOP 'NOT PASSING...'
            cont = .FALSE.
c           Check if there is already a link to this segment:
            tmp_nseg = nseg
            DO iseg1 = 1, tmp_nseg
              IF (seg(iseg1,3).EQ.1) CYCLE
c             Check both ends of the current focus segment:
              DO ilink = 1, 2
                link = .FALSE.
                DO iseg2 = 1, nseg
                  IF (iseg1.EQ.iseg2) CYCLE
                  i3 = seg(iseg1,ilink)
c                 Check both ends of the test segment:
                  i4 = seg(iseg2,1)
                  i5 = seg(iseg2,2)
                  IF ((DABS(pts(i3,1)-pts(i4,1)).LT.DTOL.AND.
     .                 DABS(pts(i3,2)-pts(i4,2)).LT.DTOL).OR.
     .                (DABS(pts(i3,1)-pts(i5,1)).LT.DTOL.AND.
     .                 DABS(pts(i3,2)-pts(i5,2)).LT.DTOL)) THEN
                    link = .TRUE.
                    IF (ilink.EQ.2) seg(iseg1,3) = 1
                    EXIT
                  ENDIF
                ENDDO
                IF (debug) WRITE(fp,*) '  PASS:',pass,iseg1,link
                IF (.NOT.link) EXIT
              ENDDO

c             Both ends of the segment are attached so start
c             looking at the next segment:
              IF (link) CYCLE

              cont = .TRUE.
              DO isrf = 1, nsurface
                IF (surface(isrf)%type.NE.VESSEL_WALL.OR.    
     .              seg(iseg1,4).EQ.isrf) CYCLE
                x1 = surface(isrf)%v(1,1)
                y1 = surface(isrf)%v(2,1)
                x2 = surface(isrf)%v(1,2)
                y2 = surface(isrf)%v(2,2)
                IF ((DABS(pts(i3,1)-x1).LT.DTOL.AND.
     .               DABS(pts(i3,2)-y1).LT.DTOL).OR.
     .              (DABS(pts(i3,1)-x2).LT.DTOL.AND.
     .               DABS(pts(i3,2)-y2).LT.DTOL)) THEN

                  WRITE(fp,*) 'Wall surface 1',isrf,
     .                   surface(isrf)%index(1:2)

                  IF (res(isrf).EQ.0.0) THEN
                    res(isrf) = opt%void_res(ivoid)
                  ELSE
                    WRITE(fp,*) 'ISRF,RES:',isrf,nsurface,res(isrf)
                    WRITE(fp,*) 'X1,Y1:',x1,y1
                    WRITE(fp,*) 'X2,Y2:',x2,y2
                    WRITE(fp,*) 'RES  :',res                  
                    STOP 'NOT READY FOR THE RES...'
                  ENDIF

                  len = DSQRT((x1 - x2)**2 + (y1 - y2)**2)
                  IF (len.GT.DBLE(res(isrf))) THEN
                    tstep = 1.0D0 / DBLE(INT(len/DBLE(res(isrf))) + 1)
                  ELSE
                    tstep = 1.0D0
                  ENDIF

                  IF (debug) THEN
                    WRITE(fp,'(A,2I6,2X,2I6,2F10.4)') 
     .                ' NEW WALL SEG:',iseg1,ilink, 
     .                surface(isrf)%index(1:2),res(isrf),tstep

                    WRITE(fp,*) '    I3    =',i3
                    WRITE(fp,*) '    PTS1,2=',pts(i3,1:2)
                    WRITE(fp,*) '    ISEG11=',pts(seg(iseg1,1),1:2)
                    WRITE(fp,*) '    ISEG12=',pts(seg(iseg1,2),1:2)
                    WRITE(fp,*) '    ISRF,RES:',isrf,res(isrf)
                    WRITE(fp,*) '    X1,Y1:',x1,y1
                    WRITE(fp,*) '    X2,Y2:',x2,y2

c                    STOP 'dfsdfsd'
                  ENDIF
              
                  DO t = 0.0D0, 0.9999999D0, tstep 
                    nseg = nseg + 1
                    seg(nseg,1) = npts + 1
                    seg(nseg,2) = npts + 2
                    seg(nseg,4) = isrf
                    seg(nseg,5) = surface(isrf)%index(1)
                    npts = npts + 1
                    pts(npts,1) = x1 + t * (x2 - x1)             ! Orientation is correct
                    pts(npts,2) = y1 + t * (y2 - y1) 
                    npts = npts + 1
                    pts(npts,1) = x1 + (t + tstep) * (x2 - x1) 
                    pts(npts,2) = y1 + (t + tstep) * (y2 - y1)
                  ENDDO
                  EXIT
                ENDIF
              ENDDO
              IF (isrf.EQ.nsurface+1) 
     .          CALL ER('ProcessVoid','No link to wall surface',*99)
            ENDDO

          ENDDO

        ELSEIF (i1.GT.0.AND.i2.GT.0) THEN
c         Select the wall segments based on the index values in I1 and I2:
          DO isrf = 1, nsurface
            IF (surface(isrf)%type.NE.VESSEL_WALL.OR.
     .          (surface(isrf)%index(1).LT.i1.OR.
     .           surface(isrf)%index(1).GT.i2)) CYCLE

            WRITE(fp,*) 'Wall surface 2',isrf,surface(isrf)%index(1)

            x1 = surface(isrf)%v(1,1)
            y1 = surface(isrf)%v(2,1)
            x2 = surface(isrf)%v(1,2)
            y2 = surface(isrf)%v(2,2)

            IF (res(isrf).EQ.0.0) THEN
              res(isrf) = opt%void_res(ivoid)
            ELSE
              WRITE(fp,*) 'ISRF,RES:',isrf,nsurface,res(isrf)
              WRITE(fp,*) 'X1,Y1:',x1,y1
              WRITE(fp,*) 'X2,Y2:',x2,y2
              WRITE(fp,*) 'RES  :',res                  
              STOP 'NOT READY FOR THE RES... (2)'
            ENDIF

            len = DSQRT((x1 - x2)**2 + (y1 - y2)**2)
            IF (len.GT.DBLE(res(isrf))) THEN
              tstep = 1.0D0 / DBLE(INT(len/DBLE(res(isrf))) + 1)
            ELSE
              tstep = 1.0D0
            ENDIF

            DO t = 0.0D0, 0.9999999D0, tstep 
              nseg = nseg + 1
              seg(nseg,1) = npts + 1
              seg(nseg,2) = npts + 2
              seg(nseg,4) = isrf
              seg(nseg,5) = surface(isrf)%index(1)
              npts = npts + 1
              pts(npts,1) = x1 + t * (x2 - x1)             ! Orientation is correct
              pts(npts,2) = y1 + t * (y2 - y1) 
              npts = npts + 1
              pts(npts,1) = x1 + (t + tstep) * (x2 - x1) 
              pts(npts,2) = y1 + (t + tstep) * (y2 - y1)
            ENDDO

          ENDDO
        ENDIF
      ENDDO  ! IVOID loop

        DO isrf = 1, nsurface
           WRITE(fp,*) 'res:',isrf,res(isrf),surface(isrf)%index(1:2)
        ENDDO

c     ------------------------------------------------------------------
c...  Check for any additional line segments to include (where 'additional'
c     is from the EIRENE nomenclature meaning extra surfaces not directly
c     related to the fluid grid): 

      DO ivoid = 1, opt%nvoid
        IF (opt%void_zone(ivoid).NE.izone) CYCLE

        i1 = opt%void_add(1,ivoid)
        i2 = opt%void_add(2,ivoid)
        IF (i1.LE.0.OR.i2.LE.0) CYCLE
        
        DO isrf = 1, nsurface
          WRITE(fp,*) 'Add check',surface(isrf)%type.NE.VESSEL_WALL,
     .                surface(isrf)%index(2),i1,i2
        
          IF ((surface(isrf)%type.NE.VESSEL_WALL.AND.
     .         surface(isrf)%type.NE.HOLE_IN_GRID).OR.
     .        (surface(isrf)%index(2).LT.i1.OR.
     .         surface(isrf)%index(2).GT.i2)) CYCLE

          x1 = surface(isrf)%v(1,1)
          y1 = surface(isrf)%v(2,1)
          x2 = surface(isrf)%v(1,2)
          y2 = surface(isrf)%v(2,2)

          IF (surface(isrf)%type.EQ.HOLE_IN_GRID) THEN
            WRITE(fp,*) 'Add hole',isrf,surface(isrf)%index(2)
            nhole = nhole + 1
            xhole(nhole) = x1
            yhole(nhole) = y1
          ELSE
            WRITE(fp,*) 'Add surface',isrf,surface(isrf)%index(2)

c        DO i1 = 1, nsurface
c           WRITE(fp,*) 'res:',i1,res(i1),surface(i1)%index(1:2)
c        ENDDO
            
            IF     (res(isrf).EQ.0.0                ) THEN
              res(isrf) = opt%void_res(ivoid)
            ELSEIF (res(isrf).NE.opt%void_res(ivoid)) THEN
              WRITE(fp,*) 'ISRF,RES:',isrf,nsurface,res(isrf)
              WRITE(fp,*) 'X1,Y1:',x1,y1
              WRITE(fp,*) 'X2,Y2:',x2,y2
              WRITE(fp,*) 'RES  :',res                  
              STOP 'NOT READY FOR THE RES... (3)'
            ENDIF
            
            len = DSQRT((x1 - x2)**2 + (y1 - y2)**2)
            IF (len.GT.DBLE(res(isrf))) THEN
              tstep = 1.0D0 / DBLE(INT(len/DBLE(res(isrf))) + 1)
            ELSE
              tstep = 1.0D0
            ENDIF
            DO t = 0.0D0, 0.9999999D0, tstep 
              nseg = nseg + 1
              seg(nseg,1) = npts + 1
              seg(nseg,2) = npts + 2
              seg(nseg,4) = isrf
              seg(nseg,5) = surface(isrf)%index(2)
              npts = npts + 1
              pts(npts,1) = x1 + t * (x2 - x1)             ! Orientation is correct
              pts(npts,2) = y1 + t * (y2 - y1) 
              npts = npts + 1
              pts(npts,1) = x1 + (t + tstep) * (x2 - x1) 
              pts(npts,2) = y1 + (t + tstep) * (y2 - y1)
            ENDDO
          ENDIF        

        ENDDO  ! ISRF loop
      ENDDO  ! IVOID loop

c     ----------------------------------------------------------------------
c...  Sort segments so that the first ??? form a coherent boundary:
      DO i2 = 1, nseg-1
        DO i3 = i2+1, nseg
          IF     (DABS(pts(seg(i2,2),1)-pts(seg(i3,1),1)).LT.DTOL.AND.
     .            DABS(pts(seg(i2,2),2)-pts(seg(i3,1),2)).LT.DTOL) THEN
c            WRITE(fp,*) 'forward'
c            WRITE(fp,*) pts(seg(i2,1),1:2)
c            WRITE(fp,*) pts(seg(i2,2),1:2)
c            WRITE(fp,*) pts(seg(i3,1),1:2)
c            WRITE(fp,*) pts(seg(i3,2),1:2)
            IF (i3.EQ.i2+1) THEN
c...          Do nothing, all okay:
              EXIT
            ELSE
              seg(0   ,1:5) = seg(i2+1,1:5)
              seg(i2+1,1:5) = seg(i3  ,1:5)
              seg(i3  ,1:5) = seg(0   ,1:5)
              EXIT
            ENDIF
          ELSEIF (DABS(pts(seg(i2,2),1)-pts(seg(i3,2),1)).LT.DTOL.AND.
     .            DABS(pts(seg(i2,2),2)-pts(seg(i3,2),2)).LT.DTOL) THEN
c            WRITE(fp,*) 'backward'
c            WRITE(fp,*) pts(seg(i2,1),1:2)
c            WRITE(fp,*) pts(seg(i2,2),1:2)
c            WRITE(fp,*) pts(seg(i3,1),1:2)
c            WRITE(fp,*) pts(seg(i3,2),1:2)
            IF (i3.EQ.i2+1) THEN
              seg(0 ,1) = seg(i3,1)
              seg(i3,1) = seg(i3,2) ! Swap the order of the points
              seg(i3,2) = seg(0 ,1)
              EXIT
            ELSE
              seg(0   ,1:5) = seg(i2+1,1:5)
              seg(i2+1,1  ) = seg(i3  ,2  )  ! Swap the order of the points
              seg(i2+1,2  ) = seg(i3  ,1  )
              seg(i2+1,3:5) = seg(i3  ,3:5)
              seg(i3  ,1:5) = seg(0   ,1:5)
              EXIT
            ENDIF
          ENDIF
        ENDDO
        IF (i3.EQ.nseg+1) 
     .    CALL ER('ProcessVoid','Zone perimeter gap detected',*99)

c       Check if the path has been closed, assuming that the first segment is a
c       valid path segment (need to keep this in mind!):
c        WRITE(fp,*) 'i2,3:',i2,i3
c        WRITE(fp,*) pts(seg(i2,1),1:2)
c        WRITE(fp,*) pts(seg(i2,2),1:2)
c        WRITE(fp,*) pts(seg(i2+1,1),1:2)
c        WRITE(fp,*) pts(seg(i2+1,2),1:2)
c        DO i4 = 1, nseg
c          WRITE(fp,*) 'segment order',i4,seg(i4,5)   
c        ENDDO

        IF (DABS(pts(seg(i2+1,2),1)-pts(seg(1,1),1)).LT.DTOL.AND.
     .      DABS(pts(seg(i2+1,2),2)-pts(seg(1,1),2)).LT.DTOL) THEN 
          WRITE(fp,*) 'CLOSURE DETECTED!',i2,i3,nseg
          EXIT
        ENDIF
      ENDDO

c...  Register whether or not all of the segments that were loaded were needed
c     to close the path, or whether there are some extra ones:    
      IF (i2.NE.nseg-1) THEN
        DO i3 = 1, nseg
          WRITE(fp,*) 'segment order',i2,seg(i3,5)   
        ENDDO
        CALL WN('ProcessVoid','Non standard wall job')
      ENDIF

c      DO i2 = 1, nseg
c        WRITE(fp,*) 'segment order',i2,seg(i2,5)   
c      ENDDO

c     IF (debug) THEN
c        DO i1 = 1, nseg
c          WRITE(fp,*) 'PERIMETER:',i1,pts(seg(i1,1),1:2)
c          WRITE(fp,*) '         :',i1,pts(seg(i1,2),1:2)
c        ENDDO
c      ENDIF


c     ------------------------------------------------------------------
c...  The perimeter of the void / zone is now defined, so look to see if 
c     any additional line segments or holes should be included:     
      tmp_nseg = i2 + 1 ! nseg
      DO ivoid = 1, opt%nvoid
        IF (opt%void_zone(ivoid).NE.izone) CYCLE

c...    Set requested triangle area:
        IF (area.EQ.0.0.AND.opt%void_res(ivoid).NE.0.0) 
     .    area = 0.5 * opt%void_res(ivoid)**2

c...    Holes:
        IF     (opt%void_hole(1,ivoid).EQ.-1.0.AND.
     .          opt%void_hole(2,ivoid).EQ.-1.0) THEN
          DO isrf = 1, nsurface
            IF (surface(isrf)%type.NE.HOLE_IN_GRID) CYCLE
            x1 = surface(isrf)%v(1,1)  
            y1 = surface(isrf)%v(2,1)  
            IF (debug) WRITE(fp,'(4X,A,3I4)')
     .        'HOLE TRYING:',isrf,surface(isrf)%index(1:2)
            IF (PointInVoid(x1,y1,tmp_nseg,seg,MAXNSEG,
     .                      npts,pts,MAXNPTS,debug,fp)) THEN
              nhole = nhole + 1
              xhole(nhole) = x1
              yhole(nhole) = y1
              IF (debug) WRITE(fp,*) '   ADDING HOLE:',nhole,x1,y1
            ENDIF
          ENDDO
        ELSEIF (opt%void_hole(1,ivoid).NE.0.0.AND.
     .          opt%void_hole(2,ivoid).NE.0.0) THEN
          nhole = nhole + 1
          xhole(nhole) = opt%void_hole(1,ivoid)
          yhole(nhole) = opt%void_hole(2,ivoid)
        ENDIF

c...    Plasma conditions for all triangles in zone -- but note that
c       this does not imply any associated recycling in EIRENE:
        IF (opt%void_ne(ivoid).GT.0.0) THEN
          ne = opt%void_ne(ivoid)
          te = opt%void_te(ivoid)
          ti = opt%void_ti(ivoid)
        ENDIF

c...    Look for wall surfaces that are inside each zone:
        IF (opt%void_add(1,ivoid).EQ.-1.AND.
     .      opt%void_add(2,ivoid).EQ.-1) THEN

          DO isrf = 1, nsurface
            IF (surface(isrf)%type    .NE.VESSEL_WALL.OR.
     .          surface(isrf)%index(2).EQ.0          .OR.
     .          res(isrf).NE.0.0) CYCLE

            x1 = surface(isrf)%v(1,1)  ! Only check the first point of each surface,
            y1 = surface(isrf)%v(2,1)  ! assuming the wall isn't ill posed...
            IF (debug) WRITE(fp,'(4X,A,3I4)')
     .        'ADD  TRYING:',isrf,surface(isrf)%index(1:2)
            IF (PointInVoid(x1,y1,tmp_nseg,seg,MAXNSEG,
     .                      npts,pts,MAXNPTS,debug,fp)) THEN
c...          The point is inside the void perimeter, so add the segment to the list:
              res(isrf) = opt%void_res(ivoid)
c...          Determine if the segment needs to be sub-divided:
              x1 = surface(isrf)%v(1,1)
              y1 = surface(isrf)%v(2,1)
              x2 = surface(isrf)%v(1,2)
              y2 = surface(isrf)%v(2,2)
              len = DSQRT((x1 - x2)**2 + (y1 - y2)**2)
              IF (len.GT.DBLE(res(isrf))) THEN
                tstep = 1.0D0 / DBLE(INT(len/DBLE(res(isrf))) + 1)
              ELSE
                tstep = 1.0D0
              ENDIF
              IF (debug) THEN
                WRITE(fp,*) ' NEW ADD SEG:',isrf
c                WRITE(fp,*) '    I3    =',i3
c                WRITE(fp,*) '    PTS1,2=',pts(i3,1:2)
c                WRITE(fp,*) '    ISEG11=',pts(seg(iseg1,1),1:2)
c                WRITE(fp,*) '    ISEG12=',pts(seg(iseg1,2),1:2)
c                WRITE(fp,*) '    ISRF,RES:',isrf,res(isrf)
c                WRITE(fp,*) '    X1,Y1:',x1,y1
c                WRITE(fp,*) '    X2,Y2:',x2,y2
c                STOP 'dfsdfsd'
              ENDIF
              DO t = 0.0D0, 0.9999999D0, tstep 
                nseg = nseg + 1
                seg(nseg,1) = npts + 1
                seg(nseg,2) = npts + 2
                seg(nseg,4) = isrf
                npts = npts + 1
                pts(npts,1) = x1 + t * (x2 - x1)             ! Orientation is correct
                pts(npts,2) = y1 + t * (y2 - y1) 
                npts = npts + 1
                pts(npts,1) = x1 + (t + tstep) * (x2 - x1) 
                pts(npts,2) = y1 + (t + tstep) * (y2 - y1)
              ENDDO
            ENDIF
          ENDDO
   
        ELSE

c          DO isrf = 1, nsurface
c            IF (surface(isrf)%type    .NE.VESSEL_WALL.OR.
c     .          surface(isrf)%index(2).LT.opt%void_add(1,ivoid).OR.
c     .          surface(isrf)%index(2).GT.opt%void_add(2,ivoid)) CYCLE
c
c            res(isrf) = opt%void_res(ivoid)
cc...        Determine if the segment needs to be sub-divided:
c            x1 = surface(isrf)%v(1,1)
c            y1 = surface(isrf)%v(2,1)
c            x2 = surface(isrf)%v(1,2)
c            y2 = surface(isrf)%v(2,2)
c            len = DSQRT((x1 - x2)**2 + (y1 - y2)**2)
c            IF (len.GT.DBLE(res(isrf))) THEN
c              tstep = 1.0D0 / DBLE(INT(len/DBLE(res(isrf))) + 1)
c            ELSE
c              tstep = 1.0D0
c            ENDIF
c            IF (debug) THEN
c              WRITE(fp,*) ' NEW FORCED ADD SEG:',isrf
cc              WRITE(fp,*) '    I3    =',i3
cc              WRITE(fp,*) '    PTS1,2=',pts(i3,1:2)
cc              WRITE(fp,*) '    ISEG11=',pts(seg(iseg1,1),1:2)
cc              WRITE(fp,*) '    ISEG12=',pts(seg(iseg1,2),1:2)
cc              WRITE(fp,*) '    ISRF,RES:',isrf,res(isrf)
cc              WRITE(fp,*) '    X1,Y1:',x1,y1
cc              WRITE(fp,*) '    X2,Y2:',x2,y2
cc              STOP 'dfsdfsd'
c            ENDIF
c            DO t = 0.0D0, 0.9999999D0, tstep 
c              nseg = nseg + 1
c              seg(nseg,1) = npts + 1
c              seg(nseg,2) = npts + 2
c              seg(nseg,4) = isrf
c              npts = npts + 1
c              pts(npts,1) = x1 + t * (x2 - x1)             ! Orientation is correct
c              pts(npts,2) = y1 + t * (y2 - y1) 
c              npts = npts + 1
c              pts(npts,1) = x1 + (t + tstep) * (x2 - x1) 
c              pts(npts,2) = y1 + (t + tstep) * (y2 - y1)
c            ENDDO
c          ENDDO

        ENDIF    
 
      ENDDO


c...  Eliminate duplicate verticies:
      DO i2 = 1, npts
        DO i3 = i2+1, npts
          IF (pts(i2,1).NE.-999.0.AND.
     .        DABS(pts(i2,1)-pts(i3,1)).LT.DTOL.AND.
     .        DABS(pts(i2,2)-pts(i3,2)).LT.DTOL) THEN
            pts(i3,1) = -999.0
            pts(i3,2) = -999.0
            DO i4 = 1, nseg
              IF (seg(i4,1).EQ.i3) seg(i4,1) = i2
              IF (seg(i4,2).EQ.i3) seg(i4,2) = i2
            ENDDO
          ENDIF
        ENDDO
      ENDDO
c...  Delete points:
      DO i2 = npts, 1, -1
        IF (pts(i2,1).EQ.-999.0) THEN
          DO i3 = i2, npts-1
            pts(i3,1) = pts(i3+1,1)
            pts(i3,2) = pts(i3+1,2)
          ENDDO
          DO i4 = 1, nseg
            IF (seg(i4,1).GE.i2) seg(i4,1) = seg(i4,1) - 1
            IF (seg(i4,2).GE.i2) seg(i4,2) = seg(i4,2) - 1
          ENDDO
          npts = npts - 1
        ENDIF
      ENDDO


      DO i1 = 1, nseg
        WRITE(fp,'(A,I6,2(2F12.6,2X))') 
     .    'VOID:',i1,pts(seg(i1,1),1),pts(seg(i1,1),2),
     .               pts(seg(i1,2),1),pts(seg(i1,2),2)
      ENDDO


      fp = 99      
      OPEN(UNIT=fp,FILE='triangle.poly',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=99)      
      WRITE(fp,*) npts,2,1,0
      DO i2 = 1, npts
c        WRITE(fp,'(I6,2F12.7)') i2,pts(i2,1),pts(i2,2)
        WRITE(fp,'(I6,2F19.14)') i2,pts(i2,1),pts(i2,2)
      ENDDO
      WRITE(fp,*) nseg,0
      DO i2 = 1, nseg
        IF (i2.EQ.nseg) THEN
          WRITE(fp,'(4I6)') i2,seg(i2,1),seg(i2,2),0
        ELSE
          WRITE(fp,'(4I6)') i2,seg(i2,1),seg(i2,2),1
        ENDIF
      ENDDO
      WRITE(fp,*) nhole
      DO i2 = 1, nhole
        WRITE(fp,'(I6,2F12.7)') i2,xhole(i2),yhole(i2)
      ENDDO
      CLOSE (fp)

c...  Call triangle:
      IF (area.EQ.0.0) area = 0.01
      WRITE(command,10) 'triangle -p -q -a',area,' -Y triangle.poly>tmp'
 10   FORMAT(A,F10.8,A)
      WRITE(eirfp,*) 'COMMAND: >'//command(1:LEN_TRIM(command))//'<'
      WRITE(0    ,*) 'COMMAND: >'//command(1:LEN_TRIM(command))//'<'
      CALL CIssue(command(1:LEN_TRIM(command)),code)
      WRITE(eirfp,*) 'RETURN_CODE:',code

      icnt = icnt + 1

      CALL ReadPolyFile_06(izone,ne,te,ti)


      RETURN
 99   WRITE(0,*) 'I3    =',i3
      WRITE(0,*) 'PTS1,2=',pts(i3,1:2)
      CALL WriteEireneTriangles
      CALL DumpGrid('PROBLEM WITH VOID')
      STOP
      END
c
c ======================================================================
c
c subroutine: ReadPolyFile
c
      SUBROUTINE ReadPolyFile_06(zone,ne,te,ti)
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER, INTENT(IN) :: zone
      REAL   , INTENT(IN) :: ne,te,ti

      INTEGER fp,i1,i2,idum1,idum2

      fp = 99      
      OPEN(UNIT=fp,FILE='triangle.1.ele',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)      
      READ(fp,*) idum1
      DO i1 = 1, idum1
        READ(fp,*,ERR=98,END=98) idum2,(tri(i1+ntri)%ver(i2),i2=1,3)
        DO i2 = 1, 3
          tri(i1+ntri)%ver(i2) = tri(i1+ntri)%ver(i2) + nver
          tri(i1+ntri)%map(i2) = -1 
          tri(i1+ntri)%sid(i2) = -1
          tri(i1+ntri)%sur(i2) =  0
        ENDDO
        tri(i1+ntri)%type = VACUUM_GRID
        tri(i1+ntri)%zone = zone
        tri(i1+ntri)%index = 0
        tri(i1+ntri)%sideindex = 0
        tri(i1+ntri)%plasma = 0.0
        tri(i1+ntri)%plasma(1) = MAX(0.0,te)
        tri(i1+ntri)%plasma(2) = MAX(0.0,ti)
        tri(i1+ntri)%plasma(3) = MAX(0.0,ne)
        tri(i1+ntri)%bfield = 0.0
        tri(i1+ntri)%efield = 0.0
      ENDDO
      CLOSE (fp) 
      ntri = ntri + idum1

      fp = 99      
      OPEN(UNIT=fp,FILE='triangle.1.node',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)      
      READ(fp,*) idum1
      DO i1 = 1, idum1
        READ(fp,*,ERR=98,END=98) idum2,ver(i1+nver,1),ver(i1+nver,2)
      ENDDO
      CLOSE (fp) 
      nver = nver + idum1


c...  Build connection map:

c...  Associate with non-default standard surface and apply reflection models
c     from former additional surfaces: 

c...  Eliminate redundant verticies:


      RETURN
 98   CALL ER('ReadPolyFile','Problems with file access',*99)
 99   STOP
      END
c
c ======================================================================
c
c subroutine: DefineEireneSurfaces
c
      INTEGER FUNCTION NewEireneSurface_06(type)
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER type,index,
     .        defiliin,defilside,defilswch,defiltor,defilcol,defilcell
      REAL    defrecyct,defrecycf

      REAL       ECH
      PARAMETER (ECH=1.602192E-19)

c...  Assign defaults to surface properties:
      defiliin  = 1
      defilside = 0
      defilswch = 0
      defiltor  = 0
      defilcell = 0
      defilcol  = 1
      defrecyct = 1.0
      defrecycf = 1.0

      nsurface = nsurface + 1
      surface(nsurface)%type   = type
      surface(nsurface)%index  = 0
      surface(nsurface)%subtype = -1
      surface(nsurface)%num    = 0
      surface(nsurface)%iliin  = defiliin
      surface(nsurface)%ilside = defilside
      surface(nsurface)%ilswch = defilswch
      surface(nsurface)%iltor  = defiltor
      surface(nsurface)%ilcell = defilcell
      surface(nsurface)%ilcol  = defilcol
      surface(nsurface)%recyct = defrecyct
      surface(nsurface)%recycf = defrecycf
      surface(nsurface)%ilspt = 0
      surface(nsurface)%isrs  = 0
      surface(nsurface)%recycs = 1.0
      surface(nsurface)%recycc = 1.0

      surface(nsurface)%hard    = 0  
      WRITE(surface(nsurface)%sector,'(256X)')  ! Not sure if this is necessary, just making sure...
      surface(nsurface)%sector = 'all'

      IF     (type.EQ.VESSEL_WALL) THEN
        surface(nsurface)%surtxt   = '* vessel wall (default)'
        surface(nsurface)%reflect = LOCAL
        surface(nsurface)%ewall = -wtemp * 1.38E-23 / ECH
        surface(nsurface)%material = wmater
c...    Assume a 2-point line segment:
        surface(nsurface)%nsur = 1
        surface(nsurface)%npts(1) = 2
        surface(nsurface)%ipts(1,1) = 1
        surface(nsurface)%ipts(2,1) = 2
        surface(nsurface)%nver = 2
      ELSEIF (type.EQ.NON_DEFAULT_STANDARD) THEN
        surface(nsurface)%surtxt  = '* non-default standard (default)'
        surface(nsurface)%reflect = GLOBAL
        surface(nsurface)%ewall    = 0.0
        surface(nsurface)%material = 0.0
        surface(nsurface)%nsur = 0
        surface(nsurface)%nver = 0
      ELSEIF (type.EQ.HOLE_IN_GRID) THEN
        surface(nsurface)%surtxt  = '* hole in triangle grid'
      ELSE
        CALL ER('NewEireneSurface_06','Invalid type',*99)
      ENDIF

      NewEireneSurface_06 = nsurface

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: RefineTriangles
c
c     *** NOT CURRENTLY IN USE ***
c
      SUBROUTINE RefineTriangles_06
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      REAL*8 CalcTriangleArea

      INTEGER i1,i2,i3,iside,zone,iside2,itri2
      LOGICAL cont
      REAL    xcen,ycen
      REAL*8  x(3),y(3),area,sidelen,maxsidelen,limit


      STOP 'RefineTriangle ROUTINE NOT IN USE'
c...  *** NOTE: VER has not been converted to REAL*8 in this routine...
      RETURN

      limit = 0.0001
      zone = 3

c...  Loop over triangles:

      cont = .TRUE.
      DO WHILE (cont)
        cont = .FALSE.

        DO i1 = ntri, 1, -1

c...      Estimate triangle center (fix this up!):
          xcen = 0.0
          ycen = 0.0
          DO i2 = 1, 3
            xcen = xcen + 0.3333 * ver(tri(i1)%ver(i2),1)
            ycen = ycen + 0.3333 * ver(tri(i1)%ver(i2),2)
          ENDDO

c...      Selection rules:
          IF (tri(i1)%type.NE.VACUUM_GRID.OR.
     .        tri(i1)%zone.NE.zone.OR.
     .        xcen.LT. 0.75.OR.
     .        ycen.GT.-0.50) CYCLE         ! *HARD*

          x(1) = DBLE(ver(tri(i1)%ver(1),1))
          y(1) = DBLE(ver(tri(i1)%ver(1),2))
          x(2) = DBLE(ver(tri(i1)%ver(2),1))
          y(2) = DBLE(ver(tri(i1)%ver(2),2))
          x(3) = DBLE(ver(tri(i1)%ver(3),1))
          y(3) = DBLE(ver(tri(i1)%ver(3),2))

          area = CalcTriangleArea(x(1),y(1),x(2),y(2),x(3),y(3))

          WRITE(eirfp,*) 'I1,AREA=',i1,ntri,area

c          STOP 'sdfsd'

          IF (area.GT.limit) THEN
c...        Split triangle in half:

            cont = .TRUE.

c...        Find best side -- the longest and where it is not neighbouring
c           a triangle of a different type or in a different zone:        
            iside = 0
            maxsidelen = 0.0
            DO i2 = 1, 3
              i3 = i2 + 1
              IF (i2.EQ.3) i3 = 1
              sidelen = DSQRT( (x(i2) - x(i3))**2 + (y(i2) - y(i3))**2)
              WRITE(eirfp,*) '.......=',i1,i2,i3


              IF (tri(i1)%map(i2).LE.0.OR. ! *TEMP* is this okay
     .            tri(tri(i1)%map(i2))%type.EQ.VACUUM_GRID.AND.
     .            tri(tri(i1)%map(i2))%zone.EQ.zone.AND.            ! *HARD*
     .            sidelen.GT.maxsidelen) THEN            
                iside = i2
                maxsidelen = sidelen
              ENDIF
            ENDDO
            IF (iside.EQ.0) 
     .        CALL ER('RefineTriangles','No can find good side',*99)

            iside2 = tri(i1)%sid(iside)
            itri2  = tri(i1)%map(iside)

            WRITE(eirfp,*) 'SPLITTING:',i1,iside
            WRITE(eirfp,*) 'SPLITTING:',itri2,iside2

c...        Split original triangle:
            i2 = iside
            i3 = i2 + 1
            IF (i2.EQ.3) i3 = 1

c...        New vertex:
            nver = nver + 1
            ver(nver,1) = SNGL(0.5 * (x(i2) + x(i3)))
            ver(nver,2) = SNGL(0.5 * (y(i2) + y(i3)))
            ver(nver,3) = 0.0
c...        New triangle:
            ntri = ntri + 1
            tri(ntri)%type = VACUUM_GRID 
            tri(ntri)%zone = tri(i1)%zone
            DO i2 = 1, 3
              tri(ntri)%map(i2) = -1   ! * NEED TO DO BETTER *
              tri(ntri)%sid(i2) = -1
              tri(ntri)%sur(i2) =  0
            ENDDO
 
            IF     (iside.EQ.1) THEN
              tri(ntri)%ver(1) = nver
              tri(ntri)%ver(2) = tri(i1)%ver(2)
              tri(ntri)%ver(3) = tri(i1)%ver(3)
c...          Modify old triangle:
              tri(i1)%ver(2) = nver              
            ELSEIF (iside.EQ.2) THEN
              tri(ntri)%ver(1) = tri(i1)%ver(1)
              tri(ntri)%ver(2) = nver
              tri(ntri)%ver(3) = tri(i1)%ver(3)
c...          Modify old triangle:
              tri(i1)%ver(3) = nver              
            ELSEIF (iside.EQ.3) THEN
              tri(ntri)%ver(1) = nver
              tri(ntri)%ver(2) = tri(i1)%ver(2)
              tri(ntri)%ver(3) = tri(i1)%ver(3)
c...          Modify old triangle:
              tri(i1)%ver(3) = nver              
            ENDIF

c...        Split neigbouring triangle:

            IF (itri2.GT.0) THEN
c...          New triangle:
              ntri = ntri + 1
              tri(ntri)%type = VACUUM_GRID
              tri(ntri)%zone = tri(i1)%zone
              DO i2 = 1, 3
                tri(ntri)%map(i2) = -1 
                tri(ntri)%sid(i2) = -1
                tri(ntri)%sur(i2) =  0
              ENDDO
 
              IF     (iside2.EQ.1) THEN
                tri(ntri)%ver(1) = nver
                tri(ntri)%ver(2) = tri(itri2)%ver(2)
                tri(ntri)%ver(3) = tri(itri2)%ver(3)
c...            Modify old triangle:
                tri(itri2)%ver(2) = nver              
              ELSEIF (iside2.EQ.2) THEN
                tri(ntri)%ver(1) = tri(itri2)%ver(1)
                tri(ntri)%ver(2) = nver
                tri(ntri)%ver(3) = tri(itri2)%ver(3)
c...            Modify old triangle:
                tri(itri2)%ver(3) = nver              
              ELSEIF (iside2.EQ.3) THEN
                tri(ntri)%ver(1) = nver
                tri(ntri)%ver(2) = tri(itri2)%ver(2)
                tri(ntri)%ver(3) = tri(itri2)%ver(3)
c...            Modify old triangle:
                tri(itri2)%ver(3) = nver              
              ENDIF
            ENDIF

c            EXIT  ! * TEMP *

          ENDIF

        ENDDO

c... *TEMP*
        CALL ProcessTriangles_06(-1)

      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: WritePolyFile
c
      SUBROUTINE WritePolyFile_06(eirntri,MAXNRS,eirtri)
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER eirntri,MAXNRS
      REAL    eirtri(MAXNRS,20)
      
      INTEGER code

      REAL*8     DTOL
c      REAL       DTOL
      PARAMETER (DTOL=1.0D-07)

      INTEGER fp,i1,i2,i3,i4,i5,v1,v2,id,ik,ir,npts,
     .        nseg,seg(0:10000,2),
     .        icnt,nhole
      LOGICAL zone,hole
      REAL    area,ne,te,ti
      REAL*8  x1,x2,y1,y2,len,t,tstep,xhole(50),yhole(50),pts(10000,2)
      character*256 command

      WRITE(eirfp,*) 'HERE IN WRITEPOLYFILE'

      icnt = 0
c...CHECK BOUNDS ON SEG AND PTS AS THEY ARE ASSIGNED

      DO i1 = 1, eirntri
        IF (eirtri(i1,1).EQ.0.0) CYCLE

        nseg = 0
        npts = 0
        nhole = 0

        ne = 0.0
        te = 0.0
        ti = 0.0

c...    Build list of line segments for passing to TRIANGLE:
        DO i2 = i1+1, eirntri 
          IF (eirtri(i2,1).NE.0.0) EXIT          

c          WRITE(eirfp,*) 'DATA:',i1,i2,eirtri(i2,2)

          IF     (eirtri(i2,2).EQ.1.0) THEN  ! PARAMETER
c...        Pull line segments from existing list of triangles, which initially 
c           will just be those from the triangularization of the standard fluid
c           code magnetic grid:

            DO i3 = 1, ntri
              DO v1 = 1, 3
                IF (tri(i3)%sur(v1).GE.NINT(eirtri(i2,3)).AND.
     .              tri(i3)%sur(v1).LE.NINT(eirtri(i2,4))) THEN 

                  nseg = nseg + 1
                  seg(nseg,1) = npts + 1
                  seg(nseg,2) = npts + 2
                  v2 = v1 + 1
                  IF (v1.EQ.3) v2 = 1
                  IF (.FALSE..AND.npts.EQ.0) THEN
                    npts = npts + 1
                    pts(npts,1) = ver(tri(i3)%ver(v2),1)  ! Switch orientation??? -- for some! 
                    pts(npts,2) = ver(tri(i3)%ver(v2),2)
                    npts = npts + 1
                    pts(npts,1) = ver(tri(i3)%ver(v1),1)
                    pts(npts,2) = ver(tri(i3)%ver(v1),2)
c                    pts(npts,1) = ver(tri(i3)%ver(v2),1)  ! Switch orientation??? -- for some! 
c                    pts(npts,2) = ver(tri(i3)%ver(v2),2)
c                    npts = npts + 1
c                    pts(npts,1) = ver(tri(i3)%ver(v1),1)
c                    pts(npts,2) = ver(tri(i3)%ver(v1),2)
                  ELSE
                    npts = npts + 1
                    pts(npts,1) = ver(tri(i3)%ver(v1),1)  ! Switch orientation??? -- for some! 
                    pts(npts,2) = ver(tri(i3)%ver(v1),2)
                    npts = npts + 1
                    pts(npts,1) = ver(tri(i3)%ver(v2),1)
                    pts(npts,2) = ver(tri(i3)%ver(v2),2)
c                    pts(npts,1) = ver(tri(i3)%ver(v1),1)  ! Switch orientation??? -- for some! 
c                    pts(npts,2) = ver(tri(i3)%ver(v1),2)
c                    npts = npts + 1
c                    pts(npts,1) = ver(tri(i3)%ver(v2),1)
c                    pts(npts,2) = ver(tri(i3)%ver(v2),2)
                  ENDIF
                ENDIF
              ENDDO              
            ENDDO

          ELSEIF (eirtri(i2,2).EQ.2.0.OR.eirtri(i2,2).EQ.3.0) THEN ! PARMETER
c...        Pull line segments from the list of additional surfaces: 

c            WRITE(eirfp,*) 'I!:',i1

            DO i3 = 1, nsurface
              IF (surface(i3)%type.EQ.VESSEL_WALL.AND.
     .            ((eirtri(i2,2).EQ.2.0.AND.
     .              surface(i3)%index(1).GE.NINT(eirtri(i2,3)).AND.
     .              surface(i3)%index(1).LE.NINT(eirtri(i2,4))).OR.
     .             (eirtri(i2,2).EQ.3.0.AND.
     .              surface(i3)%index(2).GE.NINT(eirtri(i2,3)).AND.
     .              surface(i3)%index(2).LE.NINT(eirtri(i2,4))))) THEN

                x1 = surface(i3)%v(1,1)
                y1 = surface(i3)%v(2,1)
                x2 = surface(i3)%v(1,2)
                y2 = surface(i3)%v(2,2)

                len = DSQRT((x1 - x2)**2 + (y1 - y2)**2)

                IF (eirtri(i1,3).GT.0.0.AND.
     .              len.GT.DBLE(eirtri(i1,3))) THEN
                  tstep = 1.0D0 / DBLE(INT(len/DBLE(eirtri(i1,3))) + 1)
                ELSE
                  tstep = 1.0D0
                ENDIF

c                WRITE(eirfp,*) x1,y1
c                WRITE(eirfp,*) x2,y2
                DO t = 0.0D0, 0.9999999D0, tstep 
                  nseg = nseg + 1
                  seg(nseg,1) = npts + 1
                  seg(nseg,2) = npts + 2
                  npts = npts + 1
                  pts(npts,1) = x1 + t * (x2 - x1)             ! Orientation is correct
                  pts(npts,2) = y1 + t * (y2 - y1) 
                  npts = npts + 1
                  pts(npts,1) = x1 + (t + tstep) * (x2 - x1) 
                  pts(npts,2) = y1 + (t + tstep) * (y2 - y1)
                ENDDO
              ENDIF
            ENDDO
 
          ELSEIF (eirtri(i2,2).EQ.4.0) THEN
c...        Holes:
            nhole = nhole + 1
            xhole(nhole) = eirtri(i2,3)
            yhole(nhole) = eirtri(i2,4)
          ELSEIF (eirtri(i2,2).EQ.5.0) THEN
c...        Plasma conditions for all triangles in zone:
            ne = eirtri(i2,3)            
            te = eirtri(i2,4)            
            ti = eirtri(i2,5)            
          ELSE
            CALL ER('WritePolyFile','Invalid triangle grid segment',*99)
          ENDIF

        ENDDO


c...    Eliminate duplicate verticies:
        IF (.TRUE.) THEN
          DO i2 = 1, npts
            DO i3 = i2+1, npts
              IF (pts(i2,1).NE.-999.0.AND.
     .            DABS(pts(i2,1)-pts(i3,1)).LT.DTOL.AND.
     .            DABS(pts(i2,2)-pts(i3,2)).LT.DTOL) THEN
                pts(i3,1) = -999.0
                pts(i3,2) = -999.0
                DO i4 = 1, nseg
                  IF (seg(i4,1).EQ.i3) seg(i4,1) = i2
                  IF (seg(i4,2).EQ.i3) seg(i4,2) = i2
                ENDDO
              ENDIF
            ENDDO
          ENDDO
c...      Delete points:
          DO i2 = npts, 1, -1
            IF (pts(i2,1).EQ.-999.0) THEN
              DO i3 = i2, npts-1
                pts(i3,1) = pts(i3+1,1)
                pts(i3,2) = pts(i3+1,2)
              ENDDO
              DO i4 = 1, nseg
                IF (seg(i4,1).GE.i2) seg(i4,1) = seg(i4,1) - 1
                IF (seg(i4,2).GE.i2) seg(i4,2) = seg(i4,2) - 1
              ENDDO
              npts = npts - 1
            ENDIF
          ENDDO
c...      Sort segments:
          IF (.FALSE.) THEN
            DO i2 = 1, nseg-1
              DO i3 = i2+1, nseg
                IF     (DABS(pts(seg(i2,2),1)-
     .                       pts(seg(i3,1),1)).LT.DTOL.AND.
     .                  DABS(pts(seg(i2,2),2)-
     .                       pts(seg(i3,1),2)).LT.DTOL) THEN
                  IF (i3.EQ.i2+1) THEN
c...                Do nothing, all okay:
                    EXIT
                  ELSE
                    seg(0   ,1) = seg(i2+1,1)
                    seg(0   ,2) = seg(i2+1,2)
                    seg(i2+1,1) = seg(i3  ,1)
                    seg(i2+1,2) = seg(i3  ,2)
                    seg(i3  ,1) = seg(0   ,1)
                    seg(i3  ,2) = seg(0   ,2)
                    EXIT
                  ENDIF
                ELSEIF (DABS(pts(seg(i2,2),1)-
     .                       pts(seg(i3,2),1)).LT.DTOL.AND.
     .                  DABS(pts(seg(i2,2),2)-
     .                       pts(seg(i3,2),2)).LT.DTOL) THEN
                  IF (i3.EQ.i2+1) THEN
                    seg(0 ,1) = seg(i3,1)
                    seg(i3,1) = seg(i3,2)
                    seg(i3,2) = seg(0 ,1)
                    EXIT
                  ELSE
                    seg(0   ,1) = seg(i2+1,1)
                    seg(0   ,2) = seg(i2+1,2)
                    seg(i2+1,1) = seg(i3  ,2)
                    seg(i2+1,2) = seg(i3  ,1)
                    seg(i3  ,1) = seg(0   ,1)
                    seg(i3  ,2) = seg(0   ,2)
                    EXIT
                  ENDIF
                ENDIF
              ENDDO
              IF (i3.EQ.nseg+1) THEN
                WRITE(0,*) 'WARNING: SEGMENT GAP DETECTED'
              ENDIF
            ENDDO
          ENDIF

        ENDIF

        fp = 99      
        OPEN(UNIT=fp,FILE='triangle.poly',ACCESS='SEQUENTIAL',
     .       STATUS='REPLACE',ERR=99)      
        WRITE(fp,*) npts,2,1,0
        DO i2 = 1, npts
c          WRITE(fp,'(I6,2F12.7)') i2,pts(i2,1),pts(i2,2)
          WRITE(fp,'(I6,2F19.14)') i2,pts(i2,1),pts(i2,2)
        ENDDO
        WRITE(fp,*) nseg,0
        DO i2 = 1, nseg
          IF (i2.EQ.nseg) THEN
            WRITE(fp,'(4I6)') i2,seg(i2,1),seg(i2,2),0
          ELSE
            WRITE(fp,'(4I6)') i2,seg(i2,1),seg(i2,2),1
          ENDIF
        ENDDO
        WRITE(fp,*) nhole
        DO i2 = 1, nhole
          WRITE(fp,'(I6,2F12.7)') i2,xhole(i2),yhole(i2)
        ENDDO
        CLOSE (fp)

c...    Call triangle:
        area = 0.01
        IF (eirtri(i1,3).NE.0.0) area = 0.5 * eirtri(i1,3)**2
        WRITE(command,10) 'triangle -p -q -a',area,' -Y triangle.poly'
 10     FORMAT(A,F10.8,A)
        WRITE(eirfp,*) 'COMMAND: >'//command(1:LEN_TRIM(command))//'<'
        WRITE(0    ,*) 'COMMAND: >'//command(1:LEN_TRIM(command))//'<'
        CALL CIssue(command(1:LEN_TRIM(command)),code)
        WRITE(eirfp,*) 'RETURN_CODE:',code

        icnt = icnt + 1
c        IF (icnt.EQ.1) STOP 'STOP: CHECK POLY FILES'

        CALL ReadPolyFile_06(NINT(eirtri(i1,2)),ne,te,ti)

      ENDDO

c      STOP 'TEST'

      WRITE(eirfp,*) 'DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: WriteTriangleFiles
c
      SUBROUTINE WriteEireneTriangles
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER fp,i1,i2,v1,ik1,ir1,it,idum1,itri
      LOGICAL found
      REAL    version,rdum1

      REAL tot_flux

      REAL, ALLOCATABLE :: tdata(:)      

      WRITE(eirfp,*) 'WRITING TRIANGLE FILES'

      version = 1.00

      fp = 99

      ALLOCATE(tdata(ntri))

      IF (photons.EQ.-1) THEN
c...    Load ionisation data from previous EIRENE call:
        CALL LoadTriangleData(7,0,13,0,tdata,'default')  
      ELSE
        tdata = -999.0
      ENDIF

c...  Dump triangles (for OUT, not EIRENE):
      OPEN(UNIT=fp,FILE='triangles.dat',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) ntri
      DO i1 = 1, ntri
        WRITE(fp,'(I6,3(2F14.10,2X))')
     .    i1,(ver(tri(i1)%ver(i2),1),ver(tri(i1)%ver(i2),2),i2=1,3)
      ENDDO
      CLOSE(fp)      

c...  Dump vertices:
      OPEN(UNIT=fp,FILE='objects.npco_char',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) nver
      DO i1 = 1, nver
        WRITE(fp,'(I6,3F16.10)') i1,ver(i1,1)*100.0D0,ver(i1,2)*100.0D0,
     .                           0.0
c        WRITE(fp,'(I6,3F12.6)') i1,ver(i1,1)*100.0,ver(i1,2)*100.0,0.0
      ENDDO
      CLOSE(fp)      

c...  Dump sides:
      OPEN(UNIT=fp,FILE='objects.elemente',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) ntri
      DO i1 = 1, ntri
        WRITE(fp,'(I6,3(4X,3I6),4X,2I6)') 
     .    i1,(tri(i1)%ver(v1),v1=1,3),
     .    (tri(i1)%sideindex(1,v1),v1=1,3),
     .    (tri(i1)%sideindex(2,v1),v1=1,3),
     .    tri(i1)%index(1),tri(i1)%index(2)
      ENDDO
      CLOSE(fp)      

c...  Dump connection map:
      OPEN(UNIT=fp,FILE='objects.neighbors',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) ntri
      DO i1 = 1, ntri
        WRITE(fp,'(I6,4X,3(3I6,4X),4X,2I6,2X,2I4,2X,3I4)') i1,
     .    (tri(i1)%map(v1),tri(i1)%sid(v1),tri(i1)%sur(v1),v1=1,3),
     .    tri(i1)%index(1),tri(i1)%index(2),
     .    tri(i1)%type,tri(i1)%zone,
     .    tri(i1)%sideindex(5,1:3)
      ENDDO
      CLOSE(fp)      

c...  Dump plasma data:

c...  Loading magnetic field data from idl/magnetics/b.pro for "field everywhere"
c     grids for Detlev:
      IF (.FALSE..AND..NOT.tetrahedrons) THEN
        OPEN(UNIT=fp ,FILE='grid.bfield',ACCESS='SEQUENTIAL',
     .       STATUS='OLD',ERR=95)      
        READ(fp,*)
        READ(fp,*)
        DO itri = 1, ntri
          READ(fp,*) (rdum1,i1=1,4),idum1,idum1,(rdum1,i1=1,3),
     .               tri(itri)%bfield(1),tri(itri)%bfield(3),
     .               tri(itri)%bfield(2)
        ENDDO
        CLOSE(fp)
      ENDIF

      OPEN(UNIT=fp ,FILE='objects.plasma',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
c...  Header:
      WRITE(fp,'(A,F4.2,A)') 
     .  '* VERSION ',version,' OF '//
     .  fluid_code(1:LEN_TRIM(fluid_code))//
     .  ' PLASMA FILE FOR TRIANGULAR GRID'
c      WRITE(fp,'(A)') '* VERSION 1.0 OSM PLASMA FILE FOR TRIANGULAR '//
c     .                'GRID'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* BULK PLASMA DATA'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A7,2A8,2X,A10,2X,3A10,2X,3A12)') 
     .  '* Index','Te','Ti','ne','vx',
     .  'vy','vz','Bx','By','Bz'
      WRITE(fp,'(A7,2A8,2X,A10,2X,3A10,2X,3A12)') 
     .  '*      ','(eV)','(eV)','(cm-3)','(cm s-1)',
     .  '(cm s-1)','(cm s-1)','(Tesla)','(Tesla)','(Tesla)'
      WRITE(fp,*) ntri
      DO i1 = 1, ntri
        WRITE(fp,'(I7,2F8.2,2X,1P,E10.2,2X,3E10.2,2X,
     .             4E12.4,0P,6X,2I4)') i1,
     .    tri(i1)%plasma(1),           ! Te (eV)
     .    tri(i1)%plasma(2),           ! Ti (eV)
     .    tri(i1)%plasma(3)*1.0E-06,   ! ne (cm-3)
     .    tri(i1)%plasma(4)*100.0,     ! vx (cm s-1)
     .    tri(i1)%plasma(5)*100.0,     ! vy 
     .    tri(i1)%plasma(6)*100.0,     ! vz
     .    tri(i1)%bfield(1),           ! Bx (Tesla)
     .    tri(i1)%bfield(2),           ! By
     .    tri(i1)%bfield(3),           ! Bz
     .    tdata(i1),                   ! Ionisation rate from previous run (w or w/o photons)...
     .    tri(i1)%index(1),tri(i1)%index(2)
      ENDDO

      tot_flux = 0.0

c...  Target data:
      WRITE(fp,'(A)') '* TARGET DATA'
      WRITE(fp,*) ntardat
      DO i1 = 1, ntri
        IF (tri(i1)%type.NE.MAGNETIC_GRID) CYCLE
        ik1 = tri(i1)%index(1)
        ir1 = tri(i1)%index(2) 
        DO v1 = 1, 3
          IF (tri(i1)%sur      (  v1).EQ.0.OR.
     .        tri(i1)%sideindex(2,v1).EQ.0) CYCLE

c...      Find corresponding target data, as set in ProcessFluidCode, based
c         on the fluid code cell/ring indices:
          found = .FALSE.
          DO it = 1, ntardat
            IF (tardat(it,2).NE.REAL(ik1).OR.
     .          tardat(it,3).NE.REAL(ir1)) CYCLE

            IF (.NOT.found) THEN
              found = .TRUE.         

c              IF (tardat(it,7).LT.0.0) sumflux1 = sumflux1 + tardat(it,7)  ! ...debugging...
c              IF (tardat(it,7).GT.0.0) sumflux2 = sumflux2 + tardat(it,7)

              WRITE(fp,'(I7,I6,1P,E14.6,0P,2F8.2,1P,2E10.2,0P,
     .                   F6.2,1P,E10.2,0P,6X,2I4)') 
c...            Target quantities:
     .          i1,v1,                  ! Triangle index and side index
     .          tardat(it,7 ),          ! ion flux to surface for species 1 (Amps)
     .          tardat(it,6 ),          ! Te (eV)
     .          tardat(it,8 ),          ! Ti (ev)
     .          tardat(it,9 )*1.0E-06,  ! ni (cm-3)
     .          tardat(it,10)*100.0,    ! v_para (cm s-1) (not read by EIRENE as yet)  
     .          tardat(it,11),          ! Mach no.        (not read)
     .          tardat(it,12)*1.0E-04,  ! jsat (A cm-2)   (not read)
     .          ik1,ir1                 ! Fluid grid indices, for debugging only

              tot_flux = tot_flux + ABS(tardat(it,7))
c              WRITE(eirfp,*) 'FLUX 2:',it,tardat(it,7),tot_flux
            ELSE
              CALL ER('WriteTriangleFiles','Target data appears to  '//
     .                'be over-specified',*99)
            ENDIF
          ENDDO

        ENDDO
      ENDDO
      CLOSE(fp)      

      WRITE(eirfp,*) 'TOT_FLUX:',tot_flux / 1.6022E-19

      OPEN(UNIT=fp,FILE='objects.efield',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
c...  Header:
      WRITE(fp,'(A)') '* VERSION 1.0 OSM E-FIELD FILE FOR TRIANGULAR '//
     .                'GRID'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* BULK PLASMA DATA (PART DEUX)'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A7,4A12)') 
     .  '* Index','Ex','Ey','Ez','E_pot'
      WRITE(fp,'(A7,4A12)')
     .  '*      ','V m-1','V m-1','V m-1','V'
      WRITE(fp,*) ntri
      DO i1 = 1, ntri
c...    Dump efield data:
        WRITE(fp,'(I7,1P,4E12.4,10X,0P,2I4)') i1,
     .    tri(i1)%efield(1),
     .    tri(i1)%efield(2),
     .    tri(i1)%efield(3),
     .    tri(i1)%plasma(20),  ! temporary
     .    tri(i1)%index(1),tri(i1)%index(2) 
      ENDDO
      CLOSE(fp)      

      WRITE(eirfp,*) 'DONE'

      IF (ALLOCATED(tdata)) DEALLOCATE(tdata)

      RETURN
95    WRITE(0,*) 'WRITETRIANGEFILES: B-FIELD FILE NOT FOUND'
      STOP
96    WRITE(0,*) 'WRITETRIANGEFILES: PROBLEMS WITH FILE ACCESS'
      STOP
99    STOP
      END

c
c ======================================================================
c
c subroutine: AssignPlasmaQuantities_06
c
      SUBROUTINE AssignPlasmaQuantities_06
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER i1,i2

      WRITE(eirfp,*) 'ASSIGNING PLASMA QUANTITIES'

c      DO i1 = 1, ntri
c        tri(i1)%plasma = 0.0
c        tri(i1)%bfield = 0.0
c        tri(i1)%efield = 0.0
c      ENDDO

      DO i1 = 1, ntri
        DO i2 = 1, ncell
          IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .        tri(i1)%index(1).EQ.cell(i2)%index(1).AND.
     .        tri(i1)%index(2).EQ.cell(i2)%index(2)) THEN     
            tri(i1)%plasma(1:6) = cell(i2)%plasma(1:6) 
            tri(i1)%bfield(1:4) = cell(i2)%bfield(1:4) 
            tri(i1)%efield(1:3) = cell(i2)%efield(1:3) 
            tri(i1)%plasma(20)  = cell(i2)%e_pot  
          ENDIF
        ENDDO
      ENDDO

      WRITE(eirfp,*) 'DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: NewTriangle
c
      SUBROUTINE NewTriangle_06(icell)
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER icell

      ntri = ntri + 1

      tri(ntri)%type = MAGNETIC_GRID 
      tri(ntri)%zone = -1
      tri(ntri)%index(1) = cell(icell)%index(1)
      tri(ntri)%index(2) = cell(icell)%index(2)
      tri(ntri)%sideindex(1,1:3) = 0
      tri(ntri)%sideindex(2,1:3) = 0
      tri(ntri)%sideindex(3,1:3) = 0
      tri(ntri)%sideindex(4,1:3) = 0
      tri(ntri)%sideindex(5,1:3) = 0

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: AssignVertex
c
      SUBROUTINE AssignVertex_06(ivert,ipoint,iside,ilist,xlist,ylist)
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER ivert,iside,ipoint,ilist(50,2)
      REAL*8  xlist(0:50,2),ylist(0:50,2)

      IF (ilist(ipoint,iside).EQ.0) THEN
        nver = nver + 1
        ver(nver,1) = xlist(ipoint,iside)
        ver(nver,2) = ylist(ipoint,iside)
        tri(ntri)%ver(ivert) = nver
c        IF (tri(ntri)%index(1).EQ.26.AND.
c     .      tri(ntri)%index(2).EQ.39) THEN
c          WRITE(eirfp,*) 'NEW POINT ADDED',ntri,nver
c        ENDIF
      ELSE
c        IF (tri(ntri)%index(1).EQ.26.AND.
c     .      tri(ntri)%index(2).EQ.39) THEN
c          WRITE(eirfp,*) 'ASSIGNING POINT',ntri,ilist(ipoint,iside)
c          WRITE(eirfp,*) '               ',ipoint,iside
c        ENDIF
        tri(ntri)%ver(ivert) = ilist(ipoint,iside)
      ENDIF

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c function: TrianglePosition
c
c      LOGICAL FUNCTION TrianglePosition
c      USE mod_eirene06_parameters
c      USE mod_eirene06
c      IMPLICIT none
c 
c      TrianglePosition = .FALSE.
c
c      RETURN
c 99   STOP
c      END
c
c ======================================================================
c
c function: MaxTriangleAngle
c
      REAL*8 FUNCTION MaxTriangleAngle(itri)
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none
 
      INTEGER, INTENT(IN) :: itri

      REAL*8 TriangleSideLength

      REAL*8, PARAMETER :: PI = 3.1415926536D0

      REAL*8 lenA,lenB,lenC,angA,angB,angC

      lenA = TriangleSideLength(itri,1)
      lenB = TriangleSideLength(itri,2)
      lenC = TriangleSideLength(itri,3)

      angA = DACOS((lenA**2 - lenB**2 - lenC**2) / 
     .             (-2.0D0 * lenB * lenC)) * 180.0D0 / PI 

      angB = DACOS((lenB**2 - lenA**2 - lenC**2) / 
     .             (-2.0D0 * lenA * lenC)) * 180.0D0 / PI 

      angC = DACOS((lenC**2 - lenA**2 - lenB**2) / 
     .             (-2.0D0 * lenA * lenB)) * 180.0D0 / PI 

      MaxTriangleAngle = MAX(angA,MAX(angB,angC))

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c function: MaxPolygonAngle
c
      REAL*8 FUNCTION MaxPolygonAngle(n,p)
      IMPLICIT none
 
      INTEGER, INTENT(IN) :: n
      REAL*8 , INTENT(IN) :: p(n,2)

      REAL*8, PARAMETER :: PI = 3.1415926536D0
 
      INTEGER i,j,k
      REAL*8  lenI,lenJ,lenK,angle

      MaxPolygonangle = -1.0D0

      DO i = 1, n
        j = i + 1
        k = i + 2
        IF (j.GT.n) j = j - n
        IF (k.GT.n) k = k - n
        lenI = DSQRT( (p(j,1) - p(k,1))**2 + (p(j,2) - p(k,2))**2)
        lenJ = DSQRT( (p(i,1) - p(k,1))**2 + (p(i,2) - p(k,2))**2)
        lenK = DSQRT( (p(i,1) - p(j,1))**2 + (p(i,2) - p(j,2))**2)
        angle = DACOS((lenI**2 - lenJ**2 - lenK**2) / 
     .                (-2.0D0 * lenJ * lenK)) * (180.0D0 / PI) 
        IF (angle.GT.MaxPolygonAngle) MaxPolygonAngle = angle
      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c function: TriangleSideLength
c
      REAL*8 FUNCTION TriangleSideLength(itri,iside)
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none
 
      INTEGER, INTENT(IN) :: itri,iside

      INTEGER i1,i2

      i1 = iside
      i2 = iside + 1
      IF (i2.EQ.4) i2 = 1

      TriangleSideLength = 
     .  DSQRT((ver(tri(itri)%ver(i2),1) - ver(tri(itri)%ver(i1),1))**2 + 
     .        (ver(tri(itri)%ver(i2),2) - ver(tri(itri)%ver(i1),2))**2)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: BuildFluidGridTriangles
c
      SUBROUTINE BuildFluidGridTriangles_06
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      REAL*8, PARAMETER :: DTOL = 1.0D-07

c      PARAMETER (DTOL=1.0D-07) 

c      REAL    GetMach,GetJsat,GetFlux 
      LOGICAL PointOnLine,PointInPolygon
      REAL*8  MaxTriangleAngle,MaxPolygonAngle

      REAL*8, PARAMETER :: angle_limit = 170.0D0

      INTEGER i1,i2,i3,i4,nlist(2),ilist(50,2),problem_triangle,
     .        it1(2,3),it2(2,3)
      LOGICAL output,status,large_angle,print_list,test_poly
      REAL*8  x(3),y(3),s,t,xlist(0:50,2),ylist(0:50,2),slist(0:50,2),
     .        p(3,2),p_cell(4,2),xcen,ycen,max_angle1,max_angle2


c      REAL*8, ALLOCATABLE :: xvertex(:),yvertex(:)

      WRITE(eirfp,*) 'BUILDING SUPER TRIANGLES'

      CALL OutputData(85,'Building super triangles')

      problem_triangle = -10  ! 14983 

c...  Process cells and build list of triangles and verticies:
      ntri = 0
      nver = 0

      DO i1 = 1, ncell

c        IF (.NOT.(cell(i1)%index(1).EQ.211.AND.
c     .            cell(i1)%index(2).EQ.19 )) CYCENLE

        nlist = 0
        ilist = 0

        DO i2 = 1, 2
          print_list = .FALSE.
c...      Build list of points on the 14 and 23 surfaces:
          IF (i2.EQ.1) THEN
            x(1) = cell(i1)%r(1)
            y(1) = cell(i1)%z(1)
            x(2) = cell(i1)%r(4) 
            y(2) = cell(i1)%z(4)
          ELSE
            x(1) = cell(i1)%r(2)
            y(1) = cell(i1)%z(2)
            x(2) = cell(i1)%r(3)
            y(2) = cell(i1)%z(3)
          ENDIF

c...
          nlist(i2) = nlist(i2) + 1
          xlist(nlist(i2),i2) = x(1)
          ylist(nlist(i2),i2) = y(1)
          slist(nlist(i2),i2) = 0.0D0
c...      Search all other cells for points that are on the 14 and 23 sides of the 
c         current focus cell:
          DO i3 = 1, ncell
            IF (i1.EQ.i3) CYCLE   ! Speed this up, perhaps with a quick distance check, or apply a zone above...

            output = .FALSE.
c            IF (i2.EQ.2.AND.
c     .          cell(i1)%index(1).EQ.1 .AND.
c     .          cell(i1)%index(2).EQ.62.AND.
c     .          (cell(i3)%index(2).EQ.63.OR.
c     .           .FALSE.)) output = .TRUE.
            IF (i2.EQ.1.AND.cell(i1)%index(1).EQ.30.AND.               
     .                      cell(i1)%index(2).EQ.119) THEN
c            IF (i2.EQ.1.AND.cell(i1)%index(1).EQ.20.AND.               
c     .                      cell(i1)%index(2).EQ.120.AND.
c     .                      cell(i3)%index(2).EQ.119) THEN
              output = .TRUE.
              print_list = .TRUE.
            ENDIF
            IF (output) WRITE(eirfp,*) '----------------'

            IF (cell(i3)%index(1).EQ.1) THEN  ! *THIS MAY NOT WORK FOR EDGE2D*
c...          Only check points 1 and 2 on the 1st cells of a given ring.  In theory, 
c             only points on side 23 of neighbouring cells need to be checked:
              IF (i2.EQ.1) THEN
                x(3) = cell(i3)%r(2)
                y(3) = cell(i3)%z(2)
              ELSE
                x(3) = cell(i3)%r(1)
                y(3) = cell(i3)%z(1)
              ENDIF
              IF (PointOnLine(x,y,s,t,7,output)) THEN
c              IF (PointOnLine(x,y,s,t,2,output)) THEN
                IF (output) WRITE(eirfp,*) 'MORE CELLS FOUND:',i2,i3
                nlist(i2) = nlist(i2) + 1
                xlist(nlist(i2),i2) = x(3)
                ylist(nlist(i2),i2) = y(3)
                slist(nlist(i2),i2) = s
              ENDIF
              IF (output) WRITE(eirfp,*) 'CHECK 1',i3,s,t
            ENDIF
c...        Check ...:
            IF (i2.EQ.1) THEN
              x(3) = cell(i3)%r(3)
              y(3) = cell(i3)%z(3)
            ELSE
              x(3) = cell(i3)%r(4)
              y(3) = cell(i3)%z(4)
            ENDIF
            IF (PointOnLine(x,y,s,t,7,output)) THEN
c            IF (PointOnLine(x,y,s,t,2,output)) THEN
               IF (output) WRITE(eirfp,*) 'MORE CELLS FOUND 2:',i2,i3
c             WRITE(eirfp,*) 'B:',i2,i3
              nlist(i2) = nlist(i2) + 1
              xlist(nlist(i2),i2) = x(3)
              ylist(nlist(i2),i2) = y(3)
              slist(nlist(i2),i2) = s
            ENDIF
            IF (output) WRITE(eirfp,*) 'CHECK 2',i3,s,t 
          ENDDO                                         
c...
          nlist(i2) = nlist(i2) + 1
          xlist(nlist(i2),i2) = x(2)
          ylist(nlist(i2),i2) = y(2)
          slist(nlist(i2),i2) = 1.0D0

          IF (print_list) THEN
            DO i3 = 1, nlist(i2)
              WRITE(eirfp,'(A,I3,3F15.9)') 
     .          'LIST A:',i2,xlist(i3,i2),ylist(i3,i2),slist(i3,i2)
            ENDDO
          ENDIF
        ENDDO


c...    Process list of points:
        DO i2 = 1, 2
c...      Eliminate duplicates:
          DO i3 = 1, nlist(i2)-1
            DO i4 = i3+1, nlist(i2)
              IF (slist(i3,i2).EQ.slist(i4,i2)) slist(i4,i2) = -999.0D0
            ENDDO
          ENDDO
          DO i3 = nlist(i2), 1, -1
            IF (slist(i3,i2).EQ.-999.0D0) THEN 
c              WRITE(eirfp,*) 'ELIMINATING TRIANGLE POINT!', ! This should be unnecessary...
c     .          cell(i1)%index(1),cell(i1)%index(2)
              DO i4 = i3, nlist(i2)-1
                xlist(i4,i2) = xlist(i4+1,i2)
                ylist(i4,i2) = ylist(i4+1,i2)
                slist(i4,i2) = slist(i4+1,i2)
              ENDDO
              nlist(i2) = nlist(i2) - 1
            ENDIF
          ENDDO
c...      Sort:
          i3 = 1
          DO WHILE (i3.LT.nlist(i2))
            status = .FALSE.
            DO i4 = i3+1, nlist(i2)
              IF (slist(i4,i2).LT.slist(i3,i2)) THEN
                WRITE(eirfp,*) 'SORTING TRIANGLE SIDE!',    ! This sorting should be unnecessary...
     .            cell(i1)%index(1),cell(i1)%index(2)
                status = .TRUE.
                xlist(0 ,i2) = xlist(i3,i2)
                ylist(0 ,i2) = ylist(i3,i2)
                slist(0 ,i2) = slist(i3,i2)
                xlist(i3,i2) = xlist(i4,i2)
                ylist(i3,i2) = ylist(i4,i2)
                slist(i3,i2) = slist(i4,i2)
                xlist(i4,i2) = xlist(0 ,i2)
                ylist(i4,i2) = ylist(0 ,i2)
                slist(i4,i2) = slist(0 ,i2)
              ENDIF
            ENDDO
            IF (.NOT.status) i3 = i3 + 1
          ENDDO

c          DO i3 = 1, nlist(i2)
c            WRITE(eirfp,'(A,I3,3F10.4)') 
c     .        'LIST B:',i2,xlist(i3,i2),ylist(i3,i2),slist(i3,i2)
c          ENDDO

        ENDDO

c...    Check if vertices in list are already in the vertex list:
        DO i2 = 1, 2
          DO i3 = 1, nver
            DO i4 = 1, nlist(i2)
              IF (DABS(ver(i3,1)-xlist(i4,i2)).LT.DTOL.AND.
     .            DABS(ver(i3,2)-ylist(i4,i2)).LT.DTOL) THEN 

c                IF (cell(i1)%index(1).EQ.26.AND.
c     .              cell(i1)%index(2).EQ.39) THEN
                 IF (ntri.EQ.problem_triangle-1.OR.
     .               ntri.EQ.problem_triangle-2) THEN
                  WRITE(eirfp,*) 'DUPLICATE VERTEX---',
     .              i2,i4,i2,'   ',ntri,i3
                  WRITE(eirfp,*) '                   ',
     .              ver(i3,1),xlist(i4,i2)
                  WRITE(eirfp,*) '                   ',
     .              ver(i3,2),ylist(i4,i2)
                ENDIF

                ilist(i4,i2)=i3
              ENDIF
            ENDDO
          ENDDO
        ENDDO           


c...    Build triangles:

        DO i2 = 1, MAX(nlist(1),nlist(2))-1

          IF (nlist(1).GT.i2.AND.nlist(2).GT.i2) THEN
            large_angle = .FALSE.

c...        Propose a triangle scheme, which is then checked to see if it
c           makes well posed triangles before actually building them.  This
c           is necessary for "extreme" quadrangles that appear when fitting
c           the grid to a conformal wall:

            DO i3 = 1, 4
              p_cell(i3,1) = cell(i1)%r(i3)
              p_cell(i3,2) = cell(i1)%z(i3)
            ENDDO

            it1 = RESHAPE((/ i2  ,2,  i2  ,1,  i2+1,1 /),(/2,3/)) 
            it2 = RESHAPE((/ i2+1,1,  i2+1,2,  i2  ,2 /),(/2,3/)) 

            DO i3 = 1, 3
              p(i3,1) = xlist(it1(1,i3),it1(2,i3))
              p(i3,2) = ylist(it1(1,i3),it1(2,i3))
            ENDDO
            xcen = SUM(p(1:3,1)) / 3.0D0
            ycen = SUM(p(1:3,2)) / 3.0D0
            max_angle1 = MaxPolygonAngle(3,p)
            test_poly = PointInPolygon(xcen,ycen,4,p_cell)

            DO i3 = 1, 3
              p(i3,1) = xlist(it2(1,i3),it2(2,i3))
              p(i3,2) = ylist(it2(1,i3),it2(2,i3))
            ENDDO
            xcen = SUM(p(1:3,1)) / 3.0D0
            ycen = SUM(p(1:3,2)) / 3.0D0
            max_angle1 = MAX(max_angle1,MaxPolygonAngle(3,p))
            test_poly = test_poly.AND.PointInPolygon(xcen,ycen,4,p_cell)

            IF (ABS(ntri+1-problem_triangle).LT.4) THEN
              WRITE(0,*) '------------------------'
              WRITE(0,*) 'FREAK:',ntri+1,ntri+2,
     .                       test_poly,max_angle1
            ENDIF

            IF (.NOT.test_poly.OR.max_angle1.GT.angle_limit) THEN
c...          If either test for triangle quality has failed, reverse the
c             orientation of the polygon and see if that is any better:
              it1(2,3) = 2
              it2(2,3) = 1
              IF (.NOT.test_poly) THEN
c...            The polygon test has failed so there's no need to check
c               if the change in orientation is an improvement or not, so
c               just proceed, fingers crossed:
                IF (ABS(ntri+1-problem_triangle).LT.4) THEN
                  WRITE(0,*) '   FREAK POWER:',ntri+1,ntri+2
                ENDIF
              ELSE
                DO i3 = 1, 3
                  p(i3,1) = xlist(it1(1,i3),it1(2,i3))
                  p(i3,2) = ylist(it1(1,i3),it1(2,i3))
                ENDDO
                xcen = SUM(p(1:3,1)) / 3.0D0
                ycen = SUM(p(1:3,2)) / 3.0D0
                max_angle2 = MaxPolygonAngle(3,p)
                test_poly  = PointInPolygon(xcen,ycen,4,p_cell)
                DO i3 = 1, 3
                  p(i3,1) = xlist(it2(1,i3),it2(2,i3))
                  p(i3,2) = ylist(it2(1,i3),it2(2,i3))
                ENDDO
                xcen = SUM(p(1:3,1)) / 3.0D0
                ycen = SUM(p(1:3,2)) / 3.0D0
                max_angle2 = MAX(max_angle2,MaxPolygonAngle(3,p))
                test_poly = test_poly.AND.
     .                      PointInPolygon(xcen,ycen,4,p_cell)

                IF (ABS(ntri+1-problem_triangle).LT.4) THEN
                  WRITE(0,*) '------------------------'
                  WRITE(0,*) '   FREAK AGAIN:',ntri+1,ntri+2,
     .                       test_poly,max_angle2
                ENDIF

                IF (.NOT.test_poly.OR.max_angle2.GT.max_angle1) THEN
c...              A worse situation, either because the largest interior 
c                 angle has gotten bigger or because one of the triangles
c                 is malformed, so go back to the original division scheme:
                  it1(2,3) = 1
                  it2(2,3) = 2
                  IF (ABS(ntri+1-problem_triangle).LT.4) THEN
                    WRITE(0,*) '   FREAK IN REVERSE'
                  ENDIF
                ENDIF
              ENDIF
            ENDIF

            CALL NewTriangle_06(i1)
            IF (it1(2,3).EQ.1) THEN
              tri(ntri)%sideindex(1,2) = cell(i1)%sideindex(1,4)
            ELSE
              tri(ntri)%sideindex(1,3) = cell(i1)%sideindex(1,2)
            ENDIF
            IF (i2.EQ.1) THEN
              tri(ntri)%sideindex(2,1) = cell(i1)%sideindex(2,1)
              tri(ntri)%sideindex(5,1) = cell(i1)%sideindex(5,1)
            ENDIF
            CALL AssignVertex_06(1,it1(1,1),it1(2,1),ilist,xlist,ylist)
            CALL AssignVertex_06(2,it1(1,2),it1(2,2),ilist,xlist,ylist)
            CALL AssignVertex_06(3,it1(1,3),it1(2,3),ilist,xlist,ylist)

            CALL NewTriangle_06(i1)
            IF (it2(2,3).EQ.2) THEN
              tri(ntri)%sideindex(1,2) = cell(i1)%sideindex(1,2)
            ELSE
              tri(ntri)%sideindex(1,3) = cell(i1)%sideindex(1,4)
            ENDIF
            IF (i2.EQ.nlist(1)-1.AND.nlist(1).EQ.nlist(2)) THEN
              tri(ntri)%sideindex(2,1) = cell(i1)%sideindex(2,3)
              tri(ntri)%sideindex(5,1) = cell(i1)%sideindex(5,3)
            ENDIF
            CALL AssignVertex_06(1,it2(1,1),it2(2,1),ilist,xlist,ylist)
            CALL AssignVertex_06(2,it2(1,2),it2(2,2),ilist,xlist,ylist)
            CALL AssignVertex_06(3,it2(1,3),it2(2,3),ilist,xlist,ylist)

          ELSEIF (nlist(1).GT.i2) THEN

            CALL NewTriangle_06(i1)
            tri(ntri)%sideindex(1,3) = cell(i1)%sideindex(1,4)
            IF (i2.EQ.nlist(1)-1) THEN
              tri(ntri)%sideindex(2,1) = cell(i1)%sideindex(2,3)
              tri(ntri)%sideindex(5,1) = cell(i1)%sideindex(5,3)
            ENDIF
            CALL AssignVertex_06(1,i2+1    ,1,ilist,xlist,ylist)
            CALL AssignVertex_06(2,nlist(2),2,ilist,xlist,ylist)
            CALL AssignVertex_06(3,i2      ,1,ilist,xlist,ylist)

            IF (ntri.EQ.problem_triangle) THEN
              WRITE(eirfp,*) 'TRIANGLE 3:',ntri,tri(ntri)%ver(1:3),
     .                                     MaxTriangleAngle(ntri)
            ENDIF

          ELSEIF (nlist(2).GT.i2) THEN

            CALL NewTriangle_06(i1)
            tri(ntri)%sideindex(1,2) = cell(i1)%sideindex(1,2)
            IF (i2.EQ.nlist(2)-1) THEN
              tri(ntri)%sideindex(2,1) = cell(i1)%sideindex(2,3)
              tri(ntri)%sideindex(5,1) = cell(i1)%sideindex(5,3)
            ENDIF
            CALL AssignVertex_06(1,nlist(1),1,ilist,xlist,ylist)
            CALL AssignVertex_06(2,i2+1    ,2,ilist,xlist,ylist)
            CALL AssignVertex_06(3,i2      ,2,ilist,xlist,ylist)

            IF (ntri.EQ.problem_triangle) THEN
              WRITE(eirfp,*) 'TRIANGLE 4:',ntri,tri(ntri)%ver(1:3),
     .                                     MaxTriangleAngle(ntri)
            ENDIF

          ELSE
            CALL ER('BuildFluidGridTriangles','Unknown situation',*99)
          ENDIF

        ENDDO

      ENDDO

c      CALL WriteEireneTriangles
c      CALL DumpGrid('PROBLEM WITH TRIANGLE MAP')
c      STOP 'CRAPPO!'

c      STOP 'DUMPING'
      WRITE(eirfp,*) 'DONE'

      RETURN
99    STOP
      END
c
c
c ======================================================================
c
c  subroutine: SaveTriangles
c
      subroutine SaveTriangles_06
      USE mod_interface
      USE mod_eirene06
      IMPLICIT none

      INTEGER fp,i1,i2
      REAL*8  cen(3)

      fp = 99
      OPEN(UNIT=fp,FILE='triangles.raw',ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='REPLACE',ERR=98)            
      WRITE(fp,ERR=98) 1.0,ntri,nver,nsurface
      WRITE(fp,ERR=98) (tri(i1),i1=1,ntri)
      WRITE(fp,ERR=98) ((ver(i1,i2),i2=1,3),i1=1,nver)
      WRITE(fp,ERR=98) (surface(i1),i1=1,nsurface)
      CLOSE (fp)

      WRITE(88,*) 'SAVING TRIANGLES:',ntri,nver,nsurface

      IF (.TRUE.) THEN
c...    Dump grid data to an external file:
        OPEN (UNIT=fp,FILE='objects.centre',ACCESS='SEQUENTIAL',
     .        STATUS='REPLACE')      
        WRITE(fp,*) 'SHOT: 000000   TIME: 0000'
        WRITE(fp,'(2A6,2A10)') 'IK','IR','R (m)','Z (m)'
        DO i1 = 1, ntri
          CALL GetObjCentre(i1,cen)
          WRITE(fp,'(2I6,2F10.6)') 
     .      tri(i1)%index(1),tri(i1)%index(2),
     .      SNGL(cen(1)),SNGL(cen(2))
        ENDDO
        CLOSE(fp)
      ENDIF

      IF (.TRUE.) THEN
c...    Dump grid data to an external CORTEX data file:
        CALL inOpenInterface('eirene.idl.triangles',ITF_WRITE)
        CALL inPutData(tri(1:ntri)%ver(1),'VERTEX_1','N/A')
        CALL inPutData(tri(1:ntri)%ver(2),'VERTEX_2','N/A')
        CALL inPutData(tri(1:ntri)%ver(3),'VERTEX_3','N/A')
        CALL inPutData(tri(1:ntri)%map(1),'MAP_1','N/A')
        CALL inPutData(tri(1:ntri)%map(2),'MAP_2','N/A')
        CALL inPutData(tri(1:ntri)%map(3),'MAP_3','N/A')
        CALL inPutData(tri(1:ntri)%sur(1),'SURFACE_1','N/A')
        CALL inPutData(tri(1:ntri)%sur(2),'SURFACE_2','N/A')
        CALL inPutData(tri(1:ntri)%sur(3),'SURFACE_3','N/A')
        CALL inPutData(ver(1:nver,1),'VERTEX_X','m')
        CALL inPutData(ver(1:nver,2),'VERTEX_Y','m')
        CALL inCloseInterface
      ENDIF
      
      RETURN
 98   CALL ER('SaveTriangles','Problems writing data file',*99)
 99   STOP
      END
c
c ======================================================================
c
c  subroutine: CopyBlock
c
      SUBROUTINE CopyBlock(fp1,fp2)
      IMPLICIT none
      INTEGER fp1,fp2
      CHARACTER buffer*200

      DO WHILE (.TRUE.) 
        CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') THEN
          CALL WriteLine(fp2,buffer)
        ELSE
          BACKSPACE fp1
          EXIT
        ENDIF
      ENDDO

c      RETURN
 97   RETURN
c 97   CALL ER('CopyBlock','Unexpected end of file',*99)
 98   CALL ER('CopyBlock','Problems reading template file',*99)
 99   STOP
      END
c
c ======================================================================
c
c subroutine: WriteEireneInputFile_06
c
      SUBROUTINE WriteEireneInputFile_06
      USE mod_eirene06
      USE mod_sol28_global
      IMPLICIT none

      INTEGER    MAXSTRAT,MAXSDATA
      PARAMETER (MAXSTRAT=10,MAXSDATA=10)
      COMMON /EIRCOM/
     .        eir_07nsrc,eir_07opt ,eir_07stra,
     .        eir_07ind1,
     .        eir_07ind2,eir_07ind3,
     .        eir_07wght           ,eir_07data
      INTEGER eir_07nsrc,eir_07opt ,eir_07stra(MAXSTRAT),
     .        eir_07ind1(MAXSTRAT),
     .        eir_07ind2(MAXSTRAT),eir_07ind3(MAXSTRAT)
      REAL    eir_07wght           ,eir_07data(MAXSTRAT,MAXSDATA)

      INTEGER   ik,ik1,ik2,ir,i1,i2,i3,fp1,fp2,in,icnt,
     .          add1,ilst(1024)
      LOGICAL   output,firstcall
      REAL      x0,y0,r,zaa,roa,fact
c      REAL      x0,y0,r,vcel(MAXASCDAT),zaa,roa,fact
      CHARACTER buffer*200,geostr*4

      DATA firstcall /.TRUE./
      SAVE
c
c     Check whether DIVIMP input option requests EIRENE data file:
c
c      IF (eirdata.NE.1) RETURN

      output = .FALSE.

      IF (output) WRITE(eirfp,*) 'WRITING EIRENE INPUT FILE 06'
c
c     Initialization:
      fp1   = 97
      fp2   = 98
c      fp1   = 80
c      fp2   = 81
c      fp1   = EIRIN
c      fp2   = EIROUT

c...  eirene.input is assumed to be a template input file ...:
      OPEN(UNIT=fp1,FILE='eirene.template',FORM='FORMATTED',
     .     ERR=95,STATUS='OLD')
      OPEN(UNIT=fp2,FILE='eirene.input',FORM='FORMATTED',
     .     ERR=96,STATUS='REPLACE')
c      OPEN(UNIT=fp1,FORM='FORMATTED',ERR=95,STATUS='OLD')
c      OPEN(UNIT=fp2,FORM='FORMATTED',ERR=95,STATUS='REPLACE')

c      fp2 = 0
c      REWIND(fp1)

c      CALL MS('WriteInputFile','Using xVESM to store wall data')

cc      IF (iflexopt(6).EQ.11) THEN
c        eirtemp1 = -ctargt * 1.38E-23 / ECH
c        eirtemp2 = -cwallt * 1.38E-23 / ECH
cc      ELSE
cc        eirtemp1 = ctargt * 1.38E-23 / ECH
cc        eirtemp2 = cwallt * 1.38E-23 / ECH
cc      ENDIF


c...  Correct simulated pressure gauge volumes for cyclindrical approximation:

c
c     This is a somewhat convoluted loop at the moment, which reads
c     through an existing EIRENE input data file that serves as a
c     template for the new data file being written.  Grid and neutral
c     wall specific data are substituted into the template, as well
c     as any EIRENE options/settings that are specifiable
c     from DIVIMP (such as EIRENE execution time):
c




10    CONTINUE

      CALL ReadLine(fp1,buffer,1,*50,*98)

20    CONTINUE

      IF (buffer(1:6).EQ.'*** 0.') THEN
c
c Need to remove the requirement that the template file have an intitial
c seciton labelled *** 0...
c
c       This section has been added to EIRENE and contains options
c       for the new EIRENE code that are related to the
c       generalization of the grid:
c
        WRITE(0,*) 'DEFUNCT - HALTING CODE'
        STOP


      ELSEIF (buffer(1:6).EQ.'*** 1.') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        fp06 = fp2

c        IF (eirphoton.EQ.2) THEN
c          CALL WriteLine(fp2,buffer)
c          CALL Transferline2(fp1,fp2,buffer,3)
c        ELSE
          CALL WriteBlock01_06(fp1,fp2)
22        CALL ReadLine(fp1,buffer,1,*97,*98)
          IF (buffer(1:3).NE.'***') GOTO 22
          BACKSPACE fp1
c        ENDIF

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 2.') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        CALL WriteBlock02_06

24      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 24
        BACKSPACE fp1

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 3a') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        CALL WriteBlock03a_06

26      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 26
        BACKSPACE fp1

      ELSEIF (buffer(1:6).EQ.'*** 3b'.OR.
     .        buffer(1:6).EQ.'*** 3B') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

c        CALL WriteLine(fp2,buffer)
c        CALL CopyBlock(fp1,fp2)
        CALL WriteBlock03b_06

25      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 25
        BACKSPACE fp1

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 4.') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

c        IF (eirphoton.GT.0) THEN
c          CALL WriteLine(fp2,buffer)
c          CALL Transferline2(fp1,fp2,buffer,98)
c        ELSE

        IF (.FALSE..AND.tetrahedrons) THEN  ! TETRAHEDRONS OFF
          CALL WriteLine(fp2,buffer)
          CALL CopyBlock(fp1,fp2)
        ELSE
          CALL WriteBlock04_06
41        CALL ReadLine(fp1,buffer,1,*97,*98)
          IF (buffer(1:3).NE.'***') GOTO 41
          BACKSPACE fp1
        ENDIF

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 5.') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        IF (.FALSE..AND.tetrahedrons) THEN    ! TETRAHEDRONS OFF
          CALL WriteLine(fp2,buffer)
          CALL CopyBlock(fp1,fp2)
        ELSE
          CALL WriteBlock05_06
39        CALL ReadLine(fp1,buffer,1,*97,*98)
          IF (buffer(1:3).NE.'***') GOTO 39
          BACKSPACE fp1
        ENDIF

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 6.') THEN 
        IF (output) WRITE(eirfp,*) buffer(1:6)

        IF (photons.GT.0) THEN
          CALL WriteLine(fp2,buffer)
          CALL Transferline2(fp1,fp2,buffer,14)
        ELSEIF (.FALSE..AND.tetrahedrons) THEN     ! TETRAHEDRONS OFF
          CALL WriteLine(fp2,buffer)
          CALL CopyBlock(fp1,fp2)
        ELSE
          CALL WriteBlock06_06(fp1,fp2)
        ENDIF

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 7.') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        IF (photons.GT.0) THEN
          CALL WriteLine(fp2,buffer)
          CALL Transferline2(fp1,fp2,buffer,186)
c        ELSEIF (beam.NE.0) THEN
        ELSE
          CALL WriteBlock07_06
c        ELSE
c          CALL WriteLine(fp2,buffer)
c          CALL CopyBlock(fp1,fp2)
        ENDIF

c...    Advance input stream to the start of the next input block:
40      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 40
        BACKSPACE fp1

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:7).EQ.'*** 10.') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        IF (photons.GT.0) THEN
          WRITE(eirfp,*) 'HERE IN PHOTON CODE!'
        ELSE
        ENDIF

        CALL WriteBlock10_06
c...    Advance input stream to the start of the next input block:
42      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 42
        BACKSPACE fp1
c        CALL WriteLine(fp2,buffer)
c        CALL CopyBlock(fp1,fp2)

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 11') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        IF (.FALSE..AND.tetrahedrons) THEN   ! TETRAHEDRONS OFF
          CALL WriteLine(fp2,buffer)
          CALL CopyBlock(fp1,fp2)
        ELSE
          CALL WriteBlock11_06
        ENDIF

29      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 29
        BACKSPACE fp1

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 12') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        CALL WriteLine(fp2,buffer)
        CALL Transferline2(fp1,fp2,buffer,1)

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 13') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)
 
        IF (opt_eir%ntime.NE.0) THEN
          CALL WriteBlock13_06
        ELSE
          CALL WriteLine(fp2,buffer)
          CALL Transferline2(fp1,fp2,buffer,1)
        ENDIF

30      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 30
        BACKSPACE fp1

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 14') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        CALL WriteLine(fp2,buffer)
        CALL CopyBlock(fp1,fp2)

c        CALL WriteLine(fp2,buffer)
c        CALL Transferline2(fp1,fp2,buffer,2)
c        CALL Transferline2(fp1,fp2,buffer,16)
        IF (output) WRITE(eirfp,*) 'DONE'

        GOTO 50
      ELSE
c
c       Input block identifier is not recognized, so
c       copy the entire section as is:
c
        CALL WriteLine(fp2,buffer)

        GOTO 10
      ENDIF

      CALL ReadLine(fp1,buffer,1,*97,*98)

      IF (buffer(1:3).NE.'***') THEN
        CALL ER('WriteInputFile','Invalid template format',*99)
      ENDIF

      GOTO 20

50    CONTINUE

      CLOSE (fp1)
      CLOSE (fp2)

c      IF (.FALSE..AND.eirneut.EQ.0) THEN
c        nvesm = 0
c        write(0,*)
c        write(0,*) ' TEMPORARY BLANKING OF NEUTRAL WALL DATA! '
c        write(0,*)
c      ENDIF

c      STOP 'WRITE INPUT FILE'
      RETURN
c
c     Error code:
c 
95    WRITE(0,*) 'FILE ERROR A'
96    WRITE(0,*) 'FILE ERROR B'
      STOP
97    CALL ER('WriteInputFile','Unexpected end of file',*99)
98    CALL ER('WriteInputFile','Problems reading template file',*99)
99    WRITE(50,*) '  Last line read: '
      WRITE(50,*) '  "',buffer,'"'
c99    WRITE(EROUT,*) '  Last line read: '
c      WRITE(EROUT,*) '  "',buffer,'"'
      STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock01
c
c
c
c
      SUBROUTINE WriteBlock01_06
      USE mod_eirene06
      USE mod_sol28_global
      IMPLICIT none

      INTEGER   ntime,ntime0,iteration
      DATA iteration /0/
      SAVE

      iteration = iteration + 1

      WRITE(fp06,'(A,I6)') 
     .  '*** 1. DATA FOR OPERATING MODE (OSM), CALL ',iteration

      ntime = 0
      ntime0 = 0

      nfile = 111
      IF (photons.EQ.2) nfile = 311
      IF (tetrahedrons) nfile = 0
      IF (opt_eir%ntime.NE.0) THEN
        ntime = 1
        ntime0 = 1
        nfile = 30000
        IF (time_iteration.EQ.0.AND.opt_eir%ntime.EQ.1) nfile = 10000
        time_iteration = time_iteration + 1
      ENDIF 

      IF     (niter.GE.1) THEN
c...    BGK or photons:
        WRITE(fp06,91) 2,0,time,nfile,0,niter,0,ntime,ntime0
        WRITE(fp06,91) 1,1,0,0,1,9,1,0,5  
c        WRITE(fp06,91) 1,1,0,0,1,9,0,0,5  
        WRITE(fp06,90) 'FFFFF FFFFF FFFFF FFFFF'
      ELSEIF (.TRUE.) THEN
c...    Standard (no BGK or photons):
        WRITE(fp06,91) 2,0,time,nfile,0,0,ntime,ntime0
c        WRITE(fp06,91) 2,0,time,nfile,0,1,0,ntime
c        WRITE(fp06,91) 1,1,0,0,1,9,1,0,5  ! NGSTAL=1
c        WRITE(fp06,91) 1,1,0,0,1,0,0,0,5  ! *** MEMORY SAVING ***
        WRITE(fp06,91) 1,1,0,0,1,9,1,0,5  
c        WRITE(fp06,91) 1,1,0,0,1,9,0,0,5  
        WRITE(fp06,90) 'TFFFF FFFFF FFFFF FFFFF'
c        WRITE(fp06,90) 'FFFFF FFFF'
      ELSE
        CALL ER('WriteBlock01_06','Trouble',*99)
      ENDIF

      RETURN
90    FORMAT(A)
91    FORMAT(3I6,I6.5,20(I6:))
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock02_06
c
c
c
c
      SUBROUTINE WriteBlock02_06
      USE mod_geometry
      USE mod_eirene06
      IMPLICIT none
 
      WRITE(fp06,90) '*** 2. DATA FOR STANDARD MESH (DIVIMP)'

      IF (tetrahedrons) THEN
        WRITE(fp06,91) 1,1,1,1
        WRITE(fp06,90) 'T'
        WRITE(fp06,90) 'FFFFF FTF'
        IF (nvtxmap.NE.0) THEN
          WRITE(fp06,94) nobj+1,0,0,0,0,nvtxmap
        ELSE
          WRITE(fp06,94) nobj+1,0,0,0,0,nvtx
        ENDIF
        WRITE(fp06,90) 'CASE cmod'
        WRITE(fp06,90) 'F'
        WRITE(fp06,90) 'TFFF'
        WRITE(fp06,91) 1,0
        WRITE(fp06,92) 0.0,0.0,0.0
        WRITE(fp06,90) 'F'
        WRITE(fp06,90) 'TFFF'
        WRITE(fp06,91) 1,1
        WRITE(fp06,92) 0.0,0.0,1.0,1.0,1.0
        WRITE(fp06,90) 'F'
        WRITE(fp06,91) 0
        WRITE(fp06,90) 'F'
        WRITE(fp06,91) 0
      ELSEIF (.TRUE.) THEN
        WRITE(fp06,91) 1,1,1,1
        WRITE(fp06,90) 'T'
        WRITE(fp06,90) 'FFFFF TFF'
c        WRITE(fp06,91) ntri,0,0,0,nver
        WRITE(fp06,91) ntri+1,0,0,0,nver
        WRITE(fp06,90) 'CASE cmod'
        WRITE(fp06,90) 'F'
        WRITE(fp06,90) 'TFFF'
        WRITE(fp06,91) 1,0
        WRITE(fp06,92) 0.0,0.0,0.0
        WRITE(fp06,90) 'F'
c        WRITE(fp06,90) 'TFFF'  ! cylindrical (MAST)
c        WRITE(fp06,91) 0,0,0
c        WRITE(fp06,92) 0.0,0.0,3.768*100.0
        WRITE(fp06,90) 'FTFF'  ! toroidal 
        WRITE(fp06,91) 1,1,100
        WRITE(fp06,92) 0.0,0.0,360.0
        WRITE(fp06,90) 'F'
        WRITE(fp06,91) 0
        IF (.TRUE..OR.beam.EQ.1) THEN
          WRITE(fp06,90) 'T'  ! Additional cell for beam particles to launch into
          WRITE(fp06,91) 1    !  (PB set it up this way...)
          WRITE(fp06,93) 1.0E+00
        ELSE
          WRITE(fp06,90) 'F'
          WRITE(fp06,91) 0
        ENDIF
      ELSE
        CALL ER('WriteBlock02_06','Trouble',*99)
      ENDIF

      RETURN
90    FORMAT(A)
91    FORMAT(20(I6))
94    FORMAT(20(I8))
92    FORMAT(1P,20(E12.4))
93    FORMAT(1P,20(E12.5))
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock04_06
c
c
c
c
      SUBROUTINE WriteBlock04_06
      USE mod_eirene06
      USE mod_sol28_global
      IMPLICIT none

      INTEGER ncra,ncrm,iscd1,iscd2,natmi,nreaci

      WRITE(fp06,90) '*** 4. DATA FOR SPECIES SPECIFICATION AND '//
     .               'ATOMIC PHYSICS MODULE (OSM)'

      IF     (.TRUE.) THEN

        nreaci = 36

        WRITE(fp06,90) '* ATOMIC REACTION CARDS  NREACI='
        WRITE(fp06,91) nreaci
        WRITE(fp06,95) '  1 AMJUEL H.4 2.1.5    EI ',0,1
        WRITE(fp06,95) '  2 CONST  H.4          EI ',0,1     ! Dummy ionisation rate
        WRITE(fp06,92) -200.0,0.0,0.0,0.0,0.0,0.0
        WRITE(fp06,92)    0.0,0.0,0.0					  
        WRITE(fp06,95) '  3 AMJUEL H.102.1.5    EI ',0, 1
        WRITE(fp06,95) '  4 HYDHEL H.1 3.1.8    CX ',1, 1
        WRITE(fp06,95) '  5 HYDHEL H.3 3.1.8    CX ',1, 1
c       WRITE(fp06,95) '  6 AMJUEL H.2 2.26B0   EI ',0,56  ! For iron...
        WRITE(fp06,95) '  6 METHAN H.2 2.23      EI',0,12
        WRITE(fp06,95) '  7 METHAN H.1 3.2       CX',1,12
        WRITE(fp06,95) '  7 METHAN H.3 3.2       CX',1,12
        IF (opacity.EQ.5.AND.photons.EQ.0) THEN
          WRITE(fp06,95) '  8 H2VIBR H.4 2.1.8a   RC ',0,1
        ELSE
          WRITE(fp06,95) '  8 AMJUEL H.4 2.1.8    RC ',0,1				  
        ENDIF
        WRITE(fp06,95) '  9 HYDHEL H.2 2.2.9    EI ',0, 2,0.0,0.0,0.0
        WRITE(fp06,95) ' 10 HYDHEL H.2 2.2.5    DS ',0, 2,0.0,0.0,0.0
        WRITE(fp06,95) ' 11 HYDHEL H.2 2.2.10   DS ',0, 2,0.0,0.0,0.0
        WRITE(fp06,95) ' 13 AMJUEL H.0 0.3T     EL ',1, 2				  
        WRITE(fp06,95) ' 13 AMJUEL H.1 0.3T     EL ',1, 2				  
        WRITE(fp06,95) ' 13 AMJUEL H.3 0.3T     EL ',1, 2,0.0,0.0,0.0
        WRITE(fp06,95) ' 14 AMJUEL H.4 2.2.12   EI ',0, 2
        WRITE(fp06,95) ' 15 AMJUEL H.4 2.2.11   EI ',0, 2
        WRITE(fp06,95) ' 16 AMJUEL H.4 2.2.14   EI ',0, 2
        WRITE(fp06,95) ' 17 AMJUEL H.8 2.2.14   EI ',0, 2
        WRITE(fp06,95) ' 18 AMJUEL H.3 3.2.3    CX ',1, 2
        IF (.TRUE..OR.photons.NE.0) THEN   
c          CALL ER('WriteBlock04_06','Need to check PHOTON setup',*99)
          WRITE(fp06,96) ' 19 PHOTON H.2 HLya121.5669a RC ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 20 PHOTON H.2 HLyb102.5722a RC ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 21 PHOTON H.2 HBaa656.4667a RC ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 22 PHOTON P.1 HLya121.5669a OT ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 23 PHOTON P.1 HLyb102.5722a OT ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 24 PHOTON P.1 HBaa656.4667a OT ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 25 PHOTON H.2 HLyg97.2536a  RC ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 26 PHOTON P.1 HLyg97.2536a  OT ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 27 PHOTON H.2 HLyd94.9742a  RC ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 28 PHOTON P.1 HLyd94.9742a  OT ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 29 PHOTON H.2 HLye93.7803a  RC ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 30 PHOTON P.1 HLye93.7803a  OT ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
        ENDIF
        IF (.TRUE..OR.bgk.EQ.3) THEN
c          CALL ER('WriteBlock04_06','Need to check BGK setup',*99)
          WRITE(fp06,90) ' 31 CONST  H.2           EL  2  2'				! 17 -> 31
          WRITE(fp06,92) -2.1091E+01,0.2500E+00,0.0,0.0,0.0,0.0 
          WRITE(fp06,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0 
          WRITE(fp06,90) ' 32 CONST  H.2           EL  2  2'	                        ! 18 -> 32
          WRITE(fp06,92) -2.0589E+01,0.2500E+00,0.0,0.0,0.0,0.0
          WRITE(fp06,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0
          WRITE(fp06,90) ' 33 CONST  H.2           EL  4  4'                              ! 19 -> 33
          WRITE(fp06,92) -2.0357E+01,0.2500E+00,0.0,0.0,0.0,0.0
          WRITE(fp06,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0
        ENDIF
        IF (.TRUE..OR.beam.EQ.1) THEN 
c          CALL ER('WriteBlock04_06','Need to check BEAM setup',*99)
          WRITE(fp06,96) ' 34 AMJUEL H.1 3.1.6FJ  PI ',1,1,0.0,0.0,0.0  ! For beams...
          WRITE(fp06,95) ' 35 AMJUEL H.102.1.8     RC',0,1,1.36E+01     ! Rec e cooling...
        ENDIF
        IF (.TRUE..OR.opt_eir%ilspt.EQ.5) THEN
          WRITE(fp06,95) ' 36 AMJUEL H.2 2.26B0    EI',0,56,0.0
        ENDIF            

        natmi = 1
        IF (beam.EQ.1) natmi = natmi + 1
        IF (opt_eir%ilspt.NE.0) natmi = natmi + 1

        WRITE(fp06,90) '** 4a NEUTRAL ATOMS SPECIES CARDS: NATMI='
        WRITE(fp06,91) natmi
        ncra = 3
        ncrm = 5
        IF (bgk.EQ.3) THEN
          ncra = ncra + 2
          ncrm = ncrm + 2
        ENDIF

        WRITE(fp06,94) 1,'D(n=1)  ',2,1,1,0,1,-4,0,ncra			    
        WRITE(fp06,91) 1,115,114,0  ,30000
        WRITE(fp06,93) 2.0,0.0,0.0,0.0,1.0
        IF (photons.EQ.-1) THEN
          WRITE(fp06,91) 2,115,114,0  ,30000
          WRITE(fp06,93) 999.0,0.0,0.0,0.0,1.0
        ELSE
          WRITE(fp06,91) 2,115,114,0  ,30000
          WRITE(fp06,93) 2.0,0.0,0.0,0.0,1.0
        ENDIF
        WRITE(fp06,91) 4,114,111,114,01001				    
        WRITE(fp06,93) 0.0,0.0,0.0,0.0,1.0
        IF (bgk.EQ.3) THEN
          WRITE(fp06,91) 31,214,0,0,01001,0,111
          WRITE(fp06,93) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
          WRITE(fp06,91) 33,414,0,0,01001,0,112				    
          WRITE(fp06,93) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
        ENDIF
        IF (beam.EQ.1) THEN
          WRITE(fp06,94) 2,'D_B     ',2,1,1,0,1,-4,0,4			    
          WRITE(fp06,91) 1,115,114,0  ,30000
          WRITE(fp06,93) 2.0,0.0,0.0,0.0,1.0
          WRITE(fp06,91) 2,115,114,0  ,30000
          WRITE(fp06,93) 2.0,0.0,0.0,0.0,1.0
          WRITE(fp06,91) 4,114,111,114,01001				    
          WRITE(fp06,93) 0.0,0.0,0.0,0.0,1.0
          WRITE(fp06,91) 34,114,114,  0,01001 
          WRITE(fp06,93) 0.0,4.0,0.0
        ENDIF
        IF     (opt_eir%ilspt.EQ.2) THEN
          WRITE(fp06,94) 2,'C       ',12,6,1,0,2,2,0,2
          WRITE(fp06,91) 6,115,214,0  ,00000
          WRITE(fp06,93) 2.0,0.0,0.0,0.0,1.0
          WRITE(fp06,91) 7,114,111,214,01001
          WRITE(fp06,93) 0.0,0.0,0.0,0.0,1.0
        ELSEIF (opt_eir%ilspt.EQ.5) THEN
          WRITE(fp06,94) 2,'Fe      ',56,26,1,0,2,2,0,1
          WRITE(fp06,91) 36,115,214,0,00000
          WRITE(fp06,93) 2.0,0.0,0.0,0.0,1.0
        ELSEIF (opt_eir%ilspt.EQ.3) THEN
          WRITE(fp06,94) 2,'W       ',184,74,1,0,2,2,0,1  ! This is not the correct data!
          WRITE(fp06,91) 36,115,214,0,00000
          WRITE(fp06,93) 2.0,0.0,0.0,0.0,1.0
        ELSEIF (opt_eir%ilspt.EQ.4) THEN
          WRITE(fp06,94) 2,'Be      ',9,4,1,0,2,2,0,1     ! This is not the correct data!
          WRITE(fp06,91) 36,115,214,0,00000
          WRITE(fp06,93) 2.0,0.0,0.0,0.0,1.0
        ENDIF

        WRITE(fp06,90) '** 4b NEUTRAL MOLECULES SPECIES CARDS: NMOLI='
        WRITE(fp06,91) 1							    
        WRITE(fp06,94) 1,'D2      ',4,2,2,0,1,1,0,ncrm,0,0
        WRITE(fp06,91)  9,115,113,0,0			    
        WRITE(fp06,93) -1.5400E+01,0.0,0.0,0.0,0.0
        WRITE(fp06,91) 10,115,121,0,0
        WRITE(fp06,93) -1.0500E+01,0.0,3.0,3.0,0.0
        WRITE(fp06,91) 11,115,111,114,0
        WRITE(fp06,93) -2.5000E+01,0.0,5.0,5.0,0.0
        WRITE(fp06,91) 13,114,  0,  0,01001				    
c        WRITE(fp06,92) 0.0,0.0,0.0,0.0,1.0E-10
        WRITE(fp06,92) 0.0,0.0,0.0,0.0
        WRITE(fp06,91) 18,114,111,113,01001				    
        WRITE(fp06,92) 0.0,0.0,0.0,0.0
        IF (bgk.EQ.3) THEN
          WRITE(fp06,91) 32,314,0,0,01001,0,112
          WRITE(fp06,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
          WRITE(fp06,91) 33,514,0,0,01001,0,111
          WRITE(fp06,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
        ENDIF

        WRITE(fp06,90) '**4c TEST ION SPECIES CARDS:  NIONI ION '//
     .                 'SPECIES ARE CONSIDERED, NIONI='
        WRITE(fp06,91) 1						
        WRITE(fp06,94) 1,'D2+     ',4,2,2,1,0,-1,0,3,-1		
        WRITE(fp06,91) 14,115,111,114,0
        WRITE(fp06,92) -1.0500E+01,0.0,4.3,4.3,0.0
        WRITE(fp06,91) 15,115,124,0,0
        WRITE(fp06,92) -1.5500E+01,0.0,0.25,0.25,0.0
        WRITE(fp06,91) 16,115,121,000,30000
        WRITE(fp06,92) 16.0,0.0,10.0,0.0,0.0

        IF (photons.EQ.1.OR.photons.EQ.2) THEN

          iscd1 = 214
          iscd2 = 0
          IF (bgk.EQ.3) THEN
            iscd1 = 614
            iscd2 = 4
          ENDIF

          WRITE(fp06,90) '** 4d photons'
          WRITE(fp06,91) 6			
          WRITE(fp06,94) 1,'Ba-alpha',2,1,0,0,1,+1,0,0
          WRITE(fp06,94) 2,'Ly-alpha',2,1,0,0,1,+1,0,1
          WRITE(fp06,91) 22,iscd1,( 3+iscd2)*100+14,0,0
c          WRITE(fp06,91) 22,iscd1,314,0,0
          WRITE(fp06,92) 0.0,0.0,0.0,0.0
          WRITE(fp06,94) 3,'Ly-beta ',2,1,0,0,1,+1,0,1
          WRITE(fp06,91) 23,iscd1,( 5+iscd2)*100+14,0,0
          WRITE(fp06,92) 0.0,0.0,0.0,0.0
          WRITE(fp06,94) 4,'Ly-gamma',2,1,0,0,1,+1,0,1
          WRITE(fp06,91) 26,iscd1,( 7+iscd2)*100+14,0,0
          WRITE(fp06,92) 0.0,0.0,0.0,0.0
          WRITE(fp06,94) 5,'Ly-delta',2,1,0,0,1,+1,0,1
          WRITE(fp06,91) 28,iscd1,( 9+iscd2)*100+14,0,0
          WRITE(fp06,92) 0.0,0.0,0.0,0.0
          WRITE(fp06,94) 6,'Ly-epsi ',2,1,0,0,1,+1,0,1
          WRITE(fp06,91) 30,iscd1,(11+iscd2)*100+14,0,0
          WRITE(fp06,92) 0.0,0.0,0.0,0.0
        ELSE
          WRITE(fp06,90) '** 4d photons'
          WRITE(fp06,91) 0
        ENDIF

      ELSEIF (.FALSE.) THEN

c        WRITE(fp06,90) '* ATOMIC REACTION CARDS  NREACI='
c        WRITE(fp06,91) 22
c        WRITE(fp06,90) '  1 AMJUEL H.4 2.1.5     EI  0  1'
c        WRITE(fp06,90) '  2 AMJUEL H.102.1.5     EI  0  1'
c        WRITE(fp06,90) '  3 HYDHEL H.1 3.1.8     CX  1  1'
c        WRITE(fp06,90) '  3 HYDHEL H.3 3.1.8     CX  1  1'
c        WRITE(fp06,90) '  4 AMJUEL H.4 2.2.9     EI  0  2'
c        WRITE(fp06,90) '  5 AMJUEL H.4 2.2.5     DS  0  2'
c        WRITE(fp06,90) '  6 AMJUEL H.4 2.2.10    DS  0  2'
c        WRITE(fp06,90) '  7 AMJUEL H.4 2.2.12    DS  0  2'
c        WRITE(fp06,90) '  8 AMJUEL H.4 2.2.11    DS  0  2'
c        WRITE(fp06,90) '  9 AMJUEL H.4 2.2.14    DS  0  2'
c        WRITE(fp06,90) ' 10 AMJUEL H.8 2.2.14    DS  0  2'
c        WRITE(fp06,90) ' 12 AMJUEL H.0 0.3T     EL   1  2'
c        WRITE(fp06,90) ' 12 AMJUEL H.1 0.3T     EL   1  2'
c        WRITE(fp06,90) ' 12 AMJUEL H.3 0.3T     EL   1  2 0.00000E+00'//
c     .                 ' 0.00000E+00 0.00000E+00'
c        WRITE(fp06,90) ' 13 HYDHEL H.2 2.3.9     EI  0  4'
c        WRITE(fp06,90) ' 14 METHAN H.2 2.23      EI  0 12'
c        IF (.TRUE.) THEN
c          WRITE(fp06,90) ' 15 H2VIBR H.4 2.1.8a    RC  0  1'
c          WRITE(fp06,90) ' 16 AMJUEL H.102.1.8     RC  0  1  1.3600E 01'
c        ELSE
c          WRITE(fp06,90) ' 15 AMJUEL H.4 2.1.8     RC  0  1'
c          WRITE(fp06,90) ' 16 AMJUEL H.102.1.8     RC  0  1  1.3600E 01'
c        ENDIF
c        WRITE(fp06,90) ' 17 CONST  H.2           EL  2  2'				
c        WRITE(fp06,92) -2.1091E+01,0.2500E+00,0.0,0.0,0.0,0.0
c        WRITE(fp06,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0
c        WRITE(fp06,90) ' 18 CONST  H.2           EL  2  2'				
c        WRITE(fp06,92) -2.0589E+01,0.2500E+00,0.0,0.0,0.0,0.0
c        WRITE(fp06,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0
c        WRITE(fp06,90) ' 19 CONST  H.2           EL  4  4'
c        WRITE(fp06,92) -2.0357E+01,0.2500E+00,0.0,0.0,0.0,0.0
c        WRITE(fp06,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0
c        WRITE(fp06,90) ' 20 AMJUEL H.3 3.2.3     CX  1  2'
c        WRITE(fp06,90) ' 21 AMJUEL H.9 3.1.8     CX  1  1'
c        WRITE(fp06,90) ' 22 AMJUEL H.2 3.1.8FJ   CX  1  1'

c        WRITE(fp06,90) '*NEUTRAL ATOMS SPECIES CARDS: NATMI='
c        WRITE(fp06,91) 1							    
c        ncra = 2
c        ncrm = 5
c        IF (niter.GE.1) THEN
c          ncra = ncra + 2
c          ncrm = ncrm + 2
c        ENDIF

c        WRITE(fp06,94) 1,'D       ',2,1,1,0,1,-4,0,ncra			    
c        WRITE(fp06,91) 1,115,114,0  ,30000,000		    
c        WRITE(fp06,92) 2.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        WRITE(fp06,91) 3,114,111,114,01001				    
c        WRITE(fp06,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        IF (niter.GE.1) THEN
c          WRITE(fp06,91) 17,214,0,0,01001,0,111
c          WRITE(fp06,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c          WRITE(fp06,91) 19,414,0,0,01001,0,112				    
c          WRITE(fp06,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        ENDIF

c        WRITE(fp06,90) '** 4b NEUTRAL MOLECULES SPECIES CARDS: NMOLI='
c        WRITE(fp06,91) 1							    
c        WRITE(fp06,94) 1,'D2      ',4,2,2,0,0,2,0,ncrm
c        WRITE(fp06,91) 4,115,113,0				    
c        WRITE(fp06,92) -1.5400E+01,0.0000E+00				    
c        WRITE(fp06,91) 5,115,121,000				    
c        WRITE(fp06,92) -1.0500E+01,0.0000E+00,3.0000E+00,3.0000E+00	    
c        WRITE(fp06,91) 6,115,111,114				    
c        WRITE(fp06,92) -2.5000E+01,0.0000E+00,5.0000E+00,5.0000E+00	    
c        IF (niter.GE.1) THEN
c          WRITE(fp06,91) 18,314,0,0,01001,0,112
c          WRITE(fp06,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c          WRITE(fp06,91) 19,514,0,0,01001,0,111
c          WRITE(fp06,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        ENDIF
c        WRITE(fp06,91) 20,114,111,113,01001				    
c        WRITE(fp06,92)  0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        WRITE(fp06,91) 12,114,  0,  0,01001				    
c        WRITE(fp06,92)  0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c
c        WRITE(fp06,90) '** 4c TEST ION SPECIES CARDS:  NIONI ION '//
c     .                 'SPECIES ARE CONSIDERED, NIONI='
c        WRITE(fp06,91) 1						
c        WRITE(fp06,94) 1,'D2+     ',4,2,2,1,0,-4,0,3,-1		
c        WRITE(fp06,91) 7,115,111,114                        
c        WRITE(fp06,92) -1.0400E+01,0.0000E+00,4.3000E+00,4.3000E+00
c        WRITE(fp06,91) 8,115,124,000                        
c        WRITE(fp06,92) -1.5500E+01,0.0000E+00,0.2500E+00,0.2500E+00
c        WRITE(fp06,91) 9,115,121,000,30002                  
c        WRITE(fp06,92)  1.0000E+01,0.0000E+00,0.5000E+00,0.5000E+00

c        WRITE(fp06,90) '** 4d photons'
c        WRITE(fp06,91) 0

      ELSE
        CALL ER('WriteBlock04_06','Trouble',*99)
      ENDIF

      RETURN
90    FORMAT(A)
91    FORMAT(20(I6:))
92    FORMAT(1P,20(E12.4:))
93    FORMAT(1P,20(E12.5:))
94    FORMAT(I2,1X,A8,12(I3:))
95    FORMAT(A,2I3:,1P,3E12.5)
96    FORMAT(A,2I3:,1P,3E12.4)
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock11_06
c
c
c
c
      SUBROUTINE WriteBlock10_06
      USE mod_eirene06
      USE mod_sol28_global
      IMPLICIT none

      INTEGER i

      WRITE(fp06,90) '*** 10. DATA FOR ADDITIONAL VOLUME AND SURFACE '//
     .                       'AVERAGED TALLIES - OSM'

      IF (.TRUE.) THEN
        WRITE(fp06,91) 7,0,0,0,0,opt_eir%nadspc
        WRITE(fp06,90) '** 10A. ADDITIONAL VOLUME AVERAGED '//
     .                         'TRACKLENGTH TALLIES'
        DO i = 1, 7
          WRITE(fp06,90) '     1     1     0     1'
          WRITE(fp06,90) 'dalpha'
          WRITE(fp06,90) 'dalpha                  none'
        ENDDO
        WRITE(fp06,90) '** 10B. ADDITIONAL VOLUME AVERAGED '//
     .                         'COLLISION TALLIES'
        WRITE(fp06,90) '** 10C. ALGEBRAIC TALLIES'
        WRITE(fp06,90) '** 10D.'
        WRITE(fp06,90) '** 10E.'
        WRITE(fp06,90) '** 10F. SPECTRA'
        DO i = 1, opt_eir%nadspc
          IF (opt_eir%ispsrf_ref(i)(1:6).EQ.'eirene') THEN
            WRITE(fp06,91) opt_eir%ispsrf (i), 
     .                     opt_eir%iptyp  (i), 
     .                     opt_eir%ipsp   (i), 
     .                     opt_eir%isptyp (i), 
     .                     opt_eir%nsps   (i),
     .                     opt_eir%isrfcll(i), 
     .                     opt_eir%idirec (i)
            WRITE(fp06,93) opt_eir%spcmn  (i),
     .                     opt_eir%spcmx  (i)
            IF (opt_eir%idirec(i).NE.0) 
     .        WRITE(fp06,93) opt_eir%spcvx,
     .                       opt_eir%spcvy,
     .                       opt_eir%spcvz
          ELSE
            CALL ER('WriteBlock10_06','Invalid code referenc',*99)
          ENDIF
        ENDDO
      ELSE
        CALL ER('WriteBlock11_06','Trouble',*99)
      ENDIF

      RETURN
90    FORMAT(A)
91    FORMAT(20(I6:))
92    FORMAT(1P,20(E12.4:))
93    FORMAT(1P,20(E12.5:))
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock11_06
c
c
c
c
      SUBROUTINE WriteBlock11_06
      USE mod_eirene06
      USE mod_sol28_global
      IMPLICIT none

      WRITE(fp06,90) '*** 11. DATA FOR NUMERICAL/GRAPHICAL OUTPUT (OSM)'

      IF (.TRUE.) THEN
        WRITE(fp06,90) 'FTFFF FTFFF FFTFF TTTTT TFFFF FFFFF FFFFF'
c        WRITE(fp06,90) 'FTFFF FTFFF FFTFF TTTTT T'
c        WRITE(fp06,90) 'FTFFF TTFFF FFTFF TTTTT T'  ! Turn on particle tracking TRCGRD
c        WRITE(fp06,90) 'Ftttt ttttt tttt'
        WRITE(fp06,90) 'Ttttt ttttt tttt'
        WRITE(fp06,91) 4,MIN(1,opt_eir%nadspc)
        WRITE(fp06,91) 14,0
        WRITE(fp06,91) -2,0
        WRITE(fp06,91) -3,0
        WRITE(fp06,91) -4,0
        WRITE(fp06,91) 0
        WRITE(fp06,90) 'TTFTT FFTFT ftFFF F'
        WRITE(fp06,91) 1,ntri,1,1,1,1
        WRITE(fp06,90) 'F PEI                      1 001002'
        WRITE(fp06,90) 'F LPT                      1 003008'
        WRITE(fp06,90) 'F ENTRANCE AND COVER       1 009010'	  
        WRITE(fp06,90) 'F SOUFFLET                 2 033033 038038'
        WRITE(fp06,90) 'F VERTICAL PORT            2 039061 068068'
        WRITE(fp06,90) 'F'
        WRITE(fp06,90) 'F'
        WRITE(fp06,90) 'F'
        WRITE(fp06,92) 230.0,230.0, 80.0,0.0,-750.0
        WRITE(fp06,92)  95.0, 95.0,800.0,0.0,   0.0,750.0
        WRITE(fp06,92)  45.0, 20.0
        WRITE(fp06,91) 0,0,1,2,3,4,5,6,9,0,1
c        WRITE(fp06,91) 1,20,1,2,3,4,5,6,9,0,1
        WRITE(fp06,91) 0
      ELSE
        CALL ER('WriteBlock11_06','Trouble',*99)
      ENDIF

      RETURN
90    FORMAT(A)
91    FORMAT(20(I6:))
92    FORMAT(1P,20(E12.4:))
93    FORMAT(1P,20(E12.5:))
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock03a_06
c
c
c
c
      SUBROUTINE WriteBlock03a_06
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER i1,nstsi,instsi,def_ilcell,ilcell,n
      LOGICAL warning_reported,message
      DATA    warning_reported,message /.FALSE.,.FALSE./
      SAVE


c      WRITE(eirfp,*) 'WRITING BLOCK3a'

      def_ilcell = 0

c...  Count the number of non-default standard surfaces that have been
c     defined:
      nstsi = 0
      DO i1 = 1, nsurface
c        WRITE(eirfp,*) 'NSTSI=',nstsi 
        IF (surface(i1)%type.EQ.NON_DEFAULT_STANDARD) nstsi = nstsi + 1
      ENDDO

      WRITE(fp06,90) '*** 3a. DATA FOR NON-DEFAULT SURFACES (OSM)'
      WRITE(fp06,91) nstsi
      instsi = 0
      DO i1 = 1, nsurface
        IF (surface(i1)%type.NE.NON_DEFAULT_STANDARD) CYCLE
        instsi = instsi + 1

        ilcell = def_ilcell          
        IF (surface(i1)%ilswch.EQ.10000) ilcell = 1000  ! Special for neutral beams I believe

        n = LEN_TRIM(surface(i1)%surtxt)
        WRITE(fp06,90) surface(i1)%surtxt(1:n)  ! TRIM() doesn't work...
        WRITE(fp06,91) instsi,1,1
        WRITE(fp06,91) surface(i1)%iliin ,surface(i1)%ilside,
     .                 surface(i1)%ilswch,0             ,
     .                 surface(i1)%iltor ,surface(i1)%ilcol ,
     .                 0,ilcell,0,0

        IF (.NOT.message.AND.surface(i1)%ilspt.GT.0) THEN
          message = .TRUE.
          WRITE(0,*) '*** SPUTTERING ON IN EIRENE ***',
     .      i1,surface(i1)%ilspt,surface(i1)%iliin,
     .      surface(i1)%reflect.EQ.LOCAL
        ENDIF

        IF (.NOT.warning_reported.AND.surface(i1)%iliin.EQ.2.AND.
     .      i1.NE.3) THEN  ! Need to fix this...
          warning_reported = .TRUE.
          WRITE(0,*)
          WRITE(0,*) '--------------------------------------------'
          WRITE(0,*) '  NOTE: ILIIN=2 FOUND - SURFACE FLUXES NOT  ' 
          WRITE(0,*) '        REPORTED TO OSM (BUT ARE TO DIVIMP) ' 
          WRITE(0,*) '--------------------------------------------'
          WRITE(0,*)
        ENDIF

        IF (surface(i1)%reflect.EQ.LOCAL) THEN
          WRITE(fp06,91) 1,surface(i1)%ilspt,surface(i1)%isrs
c          WRITE(fp06,91) 1,0
          WRITE(fp06,92) surface(i1)%material,surface(i1)%ewall,
     .                   0.0,0.0,0.0,2.5
c          WRITE(fp06,92) 0.0,surface(i1)%recyct,
          WRITE(fp06,92) surface(i1)%recycf,surface(i1)%recyct,
     .                   0.0,1.0,0.5,1.0
          IF (surface(i1)%ilspt.NE.0) THEN
            WRITE(fp06,92) surface(i1)%recycs,surface(i1)%recycc
          ENDIF
        ENDIF
      ENDDO

c      WRITE(eirfp,*) 'DONE'

c      STOP 'sdfsd'

90    FORMAT(A)
91    FORMAT(20(I6:))
92    FORMAT(1P,20(E12.4:))
93    FORMAT(1P,20(E12.5:))

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock03b
c
c
c
c
      SUBROUTINE WriteBlock03b_06
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none


c      DATA material / 9642., 1206., 18474., 904./


      WRITE(fp06,80) '*** 3b. DATA FOR ADDITIONAL SURFACES (OSM)'
      IF (beam.EQ.1) THEN
        WRITE(fp06,81) 2

        WRITE(fp06,80) '* testing'
        WRITE(fp06,80) '  2.0000E 00  1.0000E 00'
        WRITE(fp06,80) '     2     2     0     0     1     2     0'//
     .                 '     0     0     0'
        WRITE(fp06,80) '  1.9400E 01  1.0000E 02 -1.0000E-05'//
     .                 '  1.9400E 01 -1.0000E 02  1.0000E-05'

        WRITE(fp06,80) '* universe, absorbing'
        WRITE(fp06,80) '  0.0000E 00  1.0000E 00'
        WRITE(fp06,80) '     2     2     0     0     0     2     0'//
     .                 '     0     0     0'
        WRITE(fp06,80) ' -5.2500E 04 -2.0000E 02  0.0000E 00'//
     .                 '  0.0000E 00  1.0000E 00  1.0000E 00'
        WRITE(fp06,80) '  0.0000E 00  0.0000E 00  0.0000E 00'//
     .                 '  0.0000E 00'

      ELSE
        WRITE(fp06,81) 1
        WRITE(fp06,80) '* testing'
        WRITE(fp06,80) '  2.0000E 00  1.0000E 00'
        WRITE(fp06,80) '     2     2     0     0     1     2     0'//
     .                 '     0     0     0'
        WRITE(fp06,80) '  1.9400E 01  1.0000E 02 -1.0000E-05'//
     .                 '  1.9400E 01 -1.0000E 02  1.0000E-05'
c        WRITE(fp06,81) 0
      ENDIF


 80   FORMAT(A)
 81   FORMAT(20(I6))
 82   FORMAT(1P,20(E12.4))
 83   FORMAT(1P,20(E12.5))       


      RETURN
      STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock05
c
c
c
c
      SUBROUTINE WriteBlock05_06
      USE mod_eirene06
      USE mod_sol28_global
      IMPLICIT none

      INTEGER   fp1,fp2,nmass,i1,ibgk,nplsi,iscd2
      CHARACTER buffer*200
      CHARACTER*8 psym(2,5)

c     nmass = NINT(crmb)
      nmass = NINT(2.0) ! Replace with assignment of mod_eirene06 variable with CRMB...


      IF (nmass.NE.1.AND.nmass.NE.2) 
     .  CALL ER('WriteBlock05','EIRENE can only run with H and D',*99)

      psym(1,1) = 'H+'
      psym(2,1) = 'D+'
      psym(1,2) = 'H(B)'
      psym(2,2) = 'D(B)'
      psym(1,3) = 'H2(B)'
      psym(2,3) = 'D2(B)'
      psym(1,4) = 'HH2(B)'
      psym(2,4) = 'DD2(B)'
      psym(1,5) = 'H2H(B)'
      psym(2,5) = 'D2D(B)'

      iscd2 = 214

      nplsi = 1
      IF (photons.EQ.1.OR.photons.EQ.2) nplsi = nplsi + 13
      IF (bgk.EQ.3) nplsi = nplsi + 4
      IF (opt_eir%ilspt.NE.0) nplsi = nplsi + 1

      WRITE(fp06,90) '*** 5. DATA FOR PLASMA-BACKGROUND (OSM) 2006'
      WRITE(fp06,90) '*BULK ION SPECIES CARDS:  NPLSI ION SPECIES '//
     .               'ARE CONSIDERED, NPLSI='
      WRITE(fp06,91) nplsi
      WRITE(fp06,94) 1,psym(nmass,1),nmass,1,1,1,1,-4,0,1      ! D+
      WRITE(fp06,91) 8,115,111,0,30000
!      WRITE(fp06,92) 0.0,0.0,0.0,0.0,1.0                      ! eirsrcmul*eirscale(11)
!      WRITE(fp06,92) 35.0,0.0,0.0,0.0,torfrac   ! added 18.07.2011 -SL
      WRITE(fp06,92) 35.0,0.0,0.0,0.0,1.0                      ! eirsrcmul*eirscale(11)
c      WRITE(fp06,92) 16.0,0.0,0.0,0.0,1.0E-15                      ! No volume recombination

      IF     (opt_eir%ilspt.EQ.2) THEN
        WRITE(fp06,94) 2,'C+      ',12 ,6 ,1,1,2,2,0,0
      ELSEIF (opt_eir%ilspt.EQ.5) THEN
        WRITE(fp06,94) 2,'Fe+     ',56 ,26,1,1,2,2,0,0
      ELSEIF (opt_eir%ilspt.EQ.3) THEN
        WRITE(fp06,94) 2,'W+      ',184,74,1,1,2,2,0,0
      ELSEIF (opt_eir%ilspt.EQ.4) THEN
        WRITE(fp06,94) 2,'Be+     ',9  ,4 ,1,1,2,2,0,0
      ENDIF

      ibgk = 0

      IF (bgk.EQ.3) THEN
c...    BGK:
        ibgk = 4
        iscd2 = 614
        WRITE(fp06,94) 2,psym(nmass,2),  nmass,1,1,0,1,-1,0,0
        WRITE(fp06,94) 3,psym(nmass,3),2*nmass,2,2,0,1,-1,0,0
        WRITE(fp06,94) 4,psym(nmass,4),  nmass,1,1,0,1,-1,0,0
        WRITE(fp06,94) 5,psym(nmass,5),2*nmass,2,2,0,1,-1,0,0
      ENDIF

      IF (photons.EQ.1.OR.photons.EQ.2) THEN
c...    Photons:

        IF (photons.EQ.1) THEN
          WRITE(fp06,94) 2+ibgk,'D_1     ',nmass,1,1,0,1,-4,0
          WRITE(fp06,94) 3+ibgk,'D_2g    ',nmass,1,1,0,1,-4,0,1
          WRITE(fp06,91) 19,0,210,iscd2,0
          WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  
        ELSE
          WRITE(fp06,94) 2+ibgk,'D_1     ',nmass,1,1,0,1,-4,0,0,0,0,0,0,
     .                   'FORT.13   ',1
          WRITE(fp06,91) 2

          WRITE(fp06,94) 3+ibgk,'D_2g    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                   'COLRAD    ',1
          WRITE(fp06,91) 19,0,210,iscd2,0
          WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  
          WRITE(fp06,95) 2,4,0,'AMJUEL H.122.1.5b   OT'
        ENDIF

        WRITE(fp06,94) 4+ibgk,'D_2c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp06,91) 19,0,210,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp06,95) 1,4,0,'AMJUEL H.122.1.8b   OT'

        WRITE(fp06,94) 5+ibgk,'D_3g    ',nmass,1,1,0,1,-4,0,1
        WRITE(fp06,91) 20,0,310,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  

        WRITE(fp06,94) 6+ibgk,'D_3c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp06,91) 20,0,310,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp06,95) 1,4,0,'AMJUEL H.122.1.8a   OT'

        WRITE(fp06,94) 7+ibgk,'D_4g    ',nmass,1,1,0,1,-4,0,1
        WRITE(fp06,91) 25,0,410,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  

        WRITE(fp06,94) 8+ibgk,'D_4c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp06,91) 25,0,410,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp06,95) 1,4,0,'AMJUEL H.122.1.8c   OT'

        WRITE(fp06,94) 9+ibgk,'D_5g    ',nmass,1,1,0,1,-4,0,1
        WRITE(fp06,91) 27,0,510,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  

        WRITE(fp06,94) 10+ibgk,'D_5c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp06,91) 27,0,510,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp06,95) 1,4,0,'AMJUEL H.122.1.8d   OT'

        WRITE(fp06,94) 11+ibgk,'D_6g    ',nmass,1,1,0,1,-4,0,1
        WRITE(fp06,91) 29,0,610,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  

        WRITE(fp06,94) 12+ibgk,'D_6c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp06,91) 29,0,610,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp06,95) 1,4,0,'AMJUEL H.122.1.8e   OT'

        WRITE(fp06,94) 13+ibgk,'Lyman_a ',nmass,1,1,0,1,-4,0,0
        WRITE(fp06,94) 14+ibgk,'Lyman_b ',nmass,1,1,0,1,-4,0,0
      ENDIF

c...
      WRITE(fp06,91) 5,-5,5,5,5
      WRITE(fp06,92) 1.0, 90.0,1.5,2.5,4.3,72.0             ! Te (bogus, not actually used in eirene)
      DO i1 = 1, nplsi
        WRITE(fp06,92) 1.0,200.0,3.0,1.0,4.3,72.0,REAL(i1)  ! Ti
      ENDDO
      DO i1 = 1, nplsi
        WRITE(fp06,92) 0.0,0.0,3.0,1.0,4.6,72.0,REAL(i1)    ! ni
      ENDDO
      DO i1 = 1, nplsi
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0, 0.0,REAL(i1)    ! vi
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0, 0.0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0, 0.0
      ENDDO
      WRITE(fp06,92) 0.1,  0.1,0.0,0.0,0.0,72.0             ! B-field


c      IF (.FALSE..AND.bgk.EQ.0.AND.photons.EQ.0) THEN
c...    Standard (no BGK or photons):
c        WRITE(fp06,91) 1
c        WRITE(fp06,91) 5,-5,5,5,5
c        WRITE(fp06,92) 1.0, 90.0,1.5,2.5,4.3,72.0
c        WRITE(fp06,92) 1.0,200.0,3.0,1.0,4.3,72.0
c        WRITE(fp06,92) 0.0,  0.0,3.0,1.0,4.6,72.0
c        WRITE(fp06,92) 0.0,  0.0,0.0,0.0,0.0, 0.0
c        WRITE(fp06,92) 0.0,  0.0,0.0,0.0,0.0, 0.0
c        WRITE(fp06,92) 0.0,  0.0,0.0,0.0,0.0, 0.0
c        WRITE(fp06.1,0.0,0.0,0.0,72.0
c      ELSE
c        CALL ER('WriteBlock05','Invalid EIR_07OPT',*99)
c      ENDIF

      RETURN
90    FORMAT(A)
91    FORMAT(20(I6))
92    FORMAT(1P,20(E12.4))
94    FORMAT(I2,1X,A8,12(I3),1X,A10,1X,I2)
95    FORMAT(3I6,1X,A)

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock06
c
c
c
c
      SUBROUTINE WriteBlock06_06(fp1,fp2)
      USE mod_eirene06
      IMPLICIT none

      INTEGER   fp1,fp2
      CHARACTER buffer*200

      WRITE(fp2,90) '*** 6. DATA FOR GENERAL REFLECTION MODEL (DIVIMP)'

      WRITE(fp2,90) 'TF'

      IF (trim_data.EQ.1) THEN
c. *HARDCODED* Need to read the DIVIMP execution directory from an enviroment
c              variable:
        WRITE(fp2,90) 'PATH  ./TRIM/'
        WRITE(fp2,90) 'D_on_W'               
        WRITE(fp2,90) 'D_on_Mo'               
        WRITE(fp2,90) 'D_on_Fe'               
        WRITE(fp2,90) 'D_on_C'               
        WRITE(fp2,90) 'D_on_Be'               
      ENDIF

      WRITE(fp2,91) 1.0
      WRITE(fp2,91) 1.0
      WRITE(fp2,91) 1.0
      WRITE(fp2,91) 1.0
      WRITE(fp2,91) ermin,50.0,0.1
c      WRITE(fp2,91) 1.0,50.0,0.1

c      WRITE(fp2,91) eirermin,50.0,0.1,eirrinteg,eireinteg

c...  Advance input stream to the start of the next input block:
40    CALL ReadLine(fp1,buffer,1,*97,*98)
      IF (buffer(1:3).NE.'***') GOTO 40
      BACKSPACE fp1

      RETURN
90    FORMAT(A)
91    FORMAT(1P,10(E12.4:),0P)
97    CALL ER('WriteBlock06','Unexpected end of file'        ,*99)
98    CALL ER('WriteBlock06','Problems reading template file',*99)
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock06
c
c
c
c
      SUBROUTINE WriteBlock07_06
      USE mod_eirene06
      IMPLICIT none

      INTEGER fp1

      INTEGER   i1
      CHARACTER buffer*200

      WRITE(fp06,90) '*** 7. DATA FOR PRIMARY SOURCES OF NEUTRALS (OSM)'
      WRITE(fp06,91) nstrata
      WRITE(fp06,91) (strata(i1)%indsrc,i1=1,nstrata)
      WRITE(fp06,92) alloc
      DO i1 = 1, nstrata
        WRITE(fp06,90) strata(i1)%txtsou(1:LEN_TRIM(strata(i1)%txtsou)) 
        WRITE(fp06,90) 'FFFFF'  
        WRITE(fp06,91) strata(i1)%npts,strata(i1)%ninitl,
     .                 strata(i1)%nemods,1
        IF (strata(i1)%flux.EQ.-999.0) THEN      ! *** HACK, FOR NOW... ***
          WRITE(fp06,92) 1.0,-1.0E-10
        ELSE
          WRITE(fp06,92) strata(i1)%flux
        ENDIF
        WRITE(fp06,90) strata(i1)%species_tag
        WRITE(fp06,91) strata(i1)%nspez
        WRITE(fp06,90) strata(i1)%distrib
        WRITE(fp06,91) 1
        WRITE(fp06,91) strata(i1)%inum,strata(i1)%indim,strata(i1)%insor
        WRITE(fp06,92) strata(i1)%sorwgt,strata(i1)%sorlim,
     .                 strata(i1)%sorind,0.0,1000.0
        WRITE(fp06,91) strata(i1)%nrsor,0,0,0,strata(i1)%nasor
        WRITE(fp06,92) strata(i1)%sorad(1:6)
        WRITE(fp06,92) strata(i1)%sorene,strata(i1)%soreni,
     .                 0.0,0.0,0.0,0.0
        WRITE(fp06,92) strata(i1)%sorcos,strata(i1)%sormax
      ENDDO

c      WRITE(eirfp,*) 'STRATA:',strata(3)%indim,strata(3)%insor


      RETURN
90    FORMAT(A)
91    FORMAT(20I6)
92    FORMAT(1P,20E12.4)
97    CALL ER('WriteBlock07','Unexpected end of file'        ,*99)
98    CALL ER('WriteBlock07','Problems reading template file',*99)
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock13
c
c
c
c
      SUBROUTINE WriteBlock13_06
      USE mod_eirene06
      IMPLICIT none

      WRITE(fp06,90) '*** 13. DATA FOR TIME DEPENDENT MODE (OSM)'
      WRITE(fp06,91) 999999
      WRITE(fp06,91) 0,1
      WRITE(fp06,92) dtimv,time0
      WRITE(fp06,90) '* 13A. DATA FOR SNAPSHOT TALLIES (OSM)'
      WRITE(fp06,91) 0

      RETURN
90    FORMAT(A)
91    FORMAT(20I6)
92    FORMAT(1P,20E12.4)
97    CALL ER('WriteBlock13','Unexpected end of file'        ,*99)
98    CALL ER('WriteBlock13','Problems reading template file',*99)
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock07
c
c
c
c
c
c ======================================================================
c
c Write EIRENE geometry file:
c
c ======================================================================
c
c subroutine: WriteGeometryFile
c
c
c
c
c ======================================================================
c
c
c ======================================================================
c
c

