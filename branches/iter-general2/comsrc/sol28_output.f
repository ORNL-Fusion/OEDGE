c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE ListTargetData(fp,title)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none      

      INTEGER  , INTENT(IN) :: fp
      CHARACTER, INTENT(IN) :: title*(*)

      INTEGER     itarget,itube,ion
      CHARACTER*2 target_tag(2) 

      ion = 1

      target_tag(LO) = 'LO'
      target_tag(HI) = 'HI'

      WRITE(fp,*)
      WRITE(fp,'(A)') 'TARGET DATA:'
      WRITE(fp,'(A)') '  '//TRIM(title)
      DO itarget = LO, HI
        WRITE(fp,*)
        WRITE(fp,'(A6,A8,A16,3A10,A6,4A10,A8,2X,A)') 
     .    'TUBE','psin','jsat','ne','ni','vi','M','pe','pi','Te','Ti',
     .    'Gamma',target_tag(itarget)
        DO itube = 1, ntube
          IF (tube(itube)%type.EQ.GRD_CORE) CYCLE
          WRITE(fp,'(I6,F8.4,1P,E16.6,3E10.2,0P,F6.2,
     .               1P,2E10.2,0P,2F10.6,F8.2)')
     .      itube,tube(itube)%psin,
     .      tube(itube)%jsat  (itarget,ion),
     .      tube(itube)%ne    (itarget),
     .      tube(itube)%ni    (itarget,ion),
     .      tube(itube)%vi    (itarget,ion),
     .      tube(itube)%machno(itarget),
     .      tube(itube)%pe    (itarget),
     .      tube(itube)%pi    (itarget,ion),
     .      tube(itube)%te    (itarget),
     .      tube(itube)%ti    (itarget,ion),
     .      tube(itube)%gamma (itarget,ion)
        ENDDO
      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
c ======================================================================
c
      SUBROUTINE GenerateOutputFiles
      IMPLICIT none

      CALL DumpData_OSM('output.end','Simulation complete')

c...  Save solution:
      CALL SaveGrid('osm.raw')

      CALL SaveFluidGridGeometry
      CALL SaveWallGeometry

      CALL User_GenerateOutputFiles

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      INTEGER FUNCTION GetTube(index,mode)
      USE mod_geometry
      USE mod_sol28_params  
      USE mod_sol28_global  
      IMPLICIT none

      INTEGER, INTENT(IN) :: index,mode

      INTEGER GetObject

      INTEGER fp,i1,ic,ic1,ic2,idim1,idim2,itube
      INTEGER last_ncell,last_ntube 
      DATA    last_ncell,last_ntube /0,0/ 
 
      SAVE

      fp = logfp

      GetTube = -1 

      IF (ntube.EQ.0) RETURN

      IF (.NOT.ALLOCATED(tube_index_map).OR.tube_modified.OR.
     .    ncell.NE.last_ncell.OR.ntube.NE.last_ntube) THEN
c...    Build the map:
        WRITE(fp,*) 'MESSAGE GetTube: Building index map'
        idim1 = MAX(ncell,nobj)
        idim2 = 2
        IF (ALLOCATED(tube_index_map)) DEALLOCATE(tube_index_map)
        ALLOCATE(tube_index_map(-1:idim1,idim2))
        tube_index_map(-1,:) = -1
        tube_index_map(0:,:) =  0
        DO ic = 1, idim1
          IF (ic.LE.nobj) i1 = GetObject(ic,IND_CELL)
          DO itube = 1, ntube
            ic1 = tube(itube)%cell_index(LO)
            ic2 = tube(itube)%cell_index(HI)
            IF (ic.LE.nobj .AND.i1.GE.ic1.AND.i1.LE.ic2) 
     .        tube_index_map(ic,1) = itube             
            IF (ic.LE.ncell.AND.ic.GE.ic1.AND.ic.LE.ic2) 
     .        tube_index_map(ic,2) = itube             
          ENDDO
        ENDDO
        last_ncell = ncell
        last_ntube = ntube
        tube_modified = .FALSE.
      ENDIF
 
      SELECTCASE (mode)
        CASE(IND_OBJECT)
          IF (index.GT.nobj) CALL ER('GetTube','Object INDEX error',*99)
          GetTube = tube_index_map(index,1)
        CASE(IND_CELL)
          IF (index.GT.ncell) CALL ER('GetTube','Cell INDEX error',*99)
          GetTube = tube_index_map(index,2)
        CASE DEFAULT
          CALL ER('GetTube','Unknown MODE',*99)
      ENDSELECT       

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      INTEGER FUNCTION GetObject(index,mode)
      USE mod_geometry
      USE mod_sol28_params  
      USE mod_sol28_global  
      IMPLICIT none

      INTEGER, INTENT(IN) :: index,mode

      INTEGER idim1,idim2,iobj,i1,fp,imap
      INTEGER last_ncell 

      DATA last_ncell / -1 /

      SAVE

      fp = logfp

      GetObject = -1 

      IF (nobj.EQ.0) RETURN

      IF (.NOT.ALLOCATED(obj_index_map).OR.obj_modified.OR.
     .    ncell.NE.last_ncell) THEN
c...    Build the map:
        WRITE(fp,*) 'MESSAGE GetObject: Building index map'
        WRITE(fp,*) .NOT.ALLOCATED(obj_index_map),obj_modified,
     .              ncell.NE.last_ncell
        idim2 = OBJ_MAXNINDEX
        idim1 = 0
        DO iobj = 1, nobj  
          IF (grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID.OR.  ! Only for fluid grid objects
     .        grp(obj(iobj)%group)%type  .NE.GRP_QUADRANGLE) CYCLE
          DO i1 = 1, idim2
            idim1 = MAX(idim1,obj(iobj)%index(i1))
          ENDDO
        ENDDO    
        IF (ALLOCATED(obj_index_map)) DEALLOCATE(obj_index_map)
        ALLOCATE(obj_index_map(idim1,idim2))
        obj_index_map = 0
        DO iobj = 1, nobj
          IF (grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID.OR.
     .        grp(obj(iobj)%group)%type  .NE.GRP_QUADRANGLE) CYCLE
          DO i1 = 1, idim2
            imap = obj(iobj)%index(i1)
            IF (imap.EQ.0) CYCLE
            IF (i1.EQ.IND_IK.OR.i1.EQ.IND_IR.OR.i1.EQ.IND_IS) CYCLE
c            WRITE(fp,*) 'GETOBJECT: ',iobj,i1,imap,obj_index_map(imap,i1)
            IF (obj_index_map(imap,i1).EQ.0) THEN
              obj_index_map(obj(iobj)%index(i1),i1) = iobj
            ELSE
c...          This simple mapping routine breaks down for 3D grids
c             that use a toroidally symmetric plasma -- need to enhance...
              CALL ER('GetObject','1:1 mapping voilation',*99)
            ENDIF
          ENDDO
        ENDDO
        last_ncell = ncell
        obj_modified = .FALSE.
      ENDIF

      GetObject = obj_index_map(index,mode)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE DumpData_OSM(fname,title)
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_sol28_targets
      USE mod_sol28_wall
      IMPLICIT none

      CHARACTER*(*) fname,title 

      INTEGER GetObject,GetTube
      REAL    GetCs2

      INTEGER   fp,it,ic,ion,i,iobj,iwall,imat,itarget
      CHARACTER tube_tag(4)*4,tag*64
      REAL      cs
      REAL*8    x(10),y(10)

      tube_tag(1) = 'SOL '
      tube_tag(2) = 'PFZ '
      tube_tag(3) = 'CORE'
      tube_tag(4) = 'N/A '

      CALL User_DumpData(fp,fname,title)

      fp = 99
      OPEN(UNIT=fp,FILE=fname(1:LEN_TRIM(fname)),ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=98)            
      
      WRITE(fp,'(A)') 'TITLE: '//TRIM(title)

      WRITE(fp,*)
      WRITE(fp,10) 'NOBJ      = ',nobj     ,'NGROUP   = ',ngrp
      WRITE(fp,10) 'NSRF      = ',nsrf     ,'NVTX     = ',nvtx
      WRITE(fp,*)
      WRITE(fp,10) 'NTUBE     = ',ntube    ,'NCELL    = ',ncell
      WRITE(fp,10) 'NION      = ',nion     ,'NFIELD   = ',nfield
      WRITE(fp,10) 'NFLUID    = ',nfluid   ,'NKINETIC = ',nkinetic
      WRITE(fp,10) 'NPIN      = ',npin     ,'NPHOTON  = ',nphoton
      WRITE(fp,*)
      WRITE(fp,10) 'NWALL     = ',nwall
      WRITE(fp,10) 'NMATERIAL = ',nmaterial
      WRITE(fp,*)
      WRITE(fp,10) 'ISEP      = ',grid%isep
      WRITE(fp,10) 'ISEP2     = ',grid%isep2
      WRITE(fp,10) 'IPFZ      = ',grid%ipfz
      WRITE(fp,*)
      WRITE(fp,10) 'NXPT      = ',grid%nxpt
      DO i = 1, grid%nxpt
        WRITE(fp,11) 'IXPT      = ',grid%ixpt(i,1:2),
     .               '(',GetTube(grid%ixpt(i,1),IND_OBJECT),
     .                   GetTube(grid%ixpt(i,2),IND_OBJECT),')',
     .               'R,ZXPT    = ',grid%rxpt(i),grid%zxpt(i)
      ENDDO
  9   FORMAT(A,I8,A)
 10   FORMAT(10(A,I8,4X))
 11   FORMAT(A,2I6,2X,A,2I6,A,8X,A,2F10.6)
 12   FORMAT(10(A,F10.6,4X))

      WRITE(fp,*)
      WRITE(fp,'(A)') 'GRID DATA:'
      WRITE(fp,*)
      WRITE(fp,12) 'R0        = ',grid%r0,'Z0        = ',grid%z0

      WRITE(fp,*)
      WRITE(fp,'(A)') 'TUBE DATA:'
      WRITE(fp,*)
      IF (ntube.GT.0) THEN
        WRITE(fp,'(2A6,A5,2A9,9A10)') 
     .    'Index','Type','#','Cell_LO','Cell_HI',
     .    'psi_n','rho(m)','L (m)',
     .    'costet_LO','costet_HI',
     .    'rp_LO'    ,'rp_HI'    ,
     .    'dds_LO'   ,'dds_HI'   
        DO it = 1, ntube
          WRITE(fp,'(I6,A6,I5,2I9,11F10.4)')
     .      it,
     .      tube_tag(tube(it)%type),
     .      tube(it)%cell_index(HI)-tube(it)%cell_index(LO)+1,
     .      tube(it)%cell_index(LO:HI),
     .      tube(it)%psin,
     .      tube(it)%rho,
     .      tube(it)%smax,
     .      tube(it)%costet(1:2),
     .      tube(it)%rp    (1:2),
     .      tube(it)%dds   (1:2)
        ENDDO
      ENDIF

      WRITE(fp,*)
      WRITE(fp,'(A)') 'TARGET DATA:'
      WRITE(fp,*)
      DO itarget = 1, ntarget
        WRITE(fp,'(I4,I4,2X,A12,I4,2X,I4,2X,1024I4)') itarget,
     .    target(itarget)%location,
     .    target(itarget)%tag,
     .    target(itarget)%position,
     .    target(itarget)%nlist,
     .    (target(itarget)%ilist(i),i=1,target(itarget)%nlist)
      ENDDO

      WRITE(fp,*)
      WRITE(fp,'(A)') 'MATERIAL DATA:'
      WRITE(fp,*)
      IF (ntube.GT.0) THEN
        WRITE(fp,'(6X,A10,3A7,A12)') 
     .    'Tag','Type','A','Z','mass(AMU)'
        DO imat = 1, nmaterial
          WRITE(fp,'(I6,A10,3I7,F12.2)')
     .      imat,
     .      TRIM(material(imat)%tag),
     .      material(imat)%type,
     .      material(imat)%A,
     .      material(imat)%Z,
     .      material(imat)%mass
        ENDDO
      ENDIF

      WRITE(fp,*)
      WRITE(fp,'(A)') 'WALL DATA:'
      WRITE(fp,*)
      WRITE(fp,'(3X,5A7,A6,A5,A9,2(2X,2A11))')
     .  'Class','Group','Index','Tube','Target','Tag','Code','T (K)',
     .  'x_v1 (m)','y_v1 (m)','x_v2 (m)','y_v2 (m)'
      DO iwall = 1, nwall
        WRITE(fp,'(I3,5I7,A6,I5,F9.1,2(2X,2F11.7))') iwall,
     .    wall(iwall)%class            ,      
     .    wall(iwall)%index(WAL_GROUP ),      
     .    wall(iwall)%index(WAL_INDEX ),      
     .    wall(iwall)%index(WAL_TUBE  ),      
     .    wall(iwall)%index(WAL_TARGET),     
     .    TRIM(wall(iwall)%material_tag),
     .    wall(iwall)%material,
     .    wall(iwall)%temperature,
     .    wall(iwall)%v1(1:2),
     .    wall(iwall)%v2(1:2)
      ENDDO

      WRITE(fp,*)
      WRITE(fp,'(A)') 'GEOMETRY DATA:'
      IF (ntube.GT.0.AND.ncell.GT.0.AND.nfield.GT.0) THEN
        DO it = 1, ntube
          DO ion = 1, nion
            WRITE(fp,*)
            WRITE(fp,10) '   TUBE =',it
            WRITE(fp,10) '   ION  =',ion
            WRITE(fp,'(3X,A8,9A11)') 
     .        'Cell','bratio','s','sbnd1','sbnd2','p','q','R','Z','vol'
            WRITE(fp,'(11X,9A11)') 
     .        ' ','(m)','(m)','(m)','(m)',' ','(m)','(m)','(m3)'
            WRITE(fp,'(11X,6F11.5,2F11.5)')
     .         tube(it)%bratio(LO),
     .         0.0,
     .        -1.0,
     .        -1.0,
     .         0.0,
     .        -1.0,
     .        -1.0,
     .        -1.0
            i = 0
            DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)        
              i = i + 1
              WRITE(fp,'(I3,I8,6F11.5,3F11.5)') i,ic,
     .          field(ic)%bratio,
     .          cell (ic)%s,
     .          cell (ic)%sbnd(1:2),
     .          cell (ic)%p,
     .          -1.0,
     .          cell (ic)%cencar(1:2),
     .          cell (ic)%vol
            ENDDO
            WRITE(fp,'(11X,6F11.5,2F11.5)')
     .        tube(it)%bratio(HI),
     .        tube(it)%smax,
     .        -1.0,
     .        -1.0,
     .        tube(it)%pmax,
     .        -1.0,
     .        -1.0,
     .        -1.0
          ENDDO
        ENDDO
      ENDIF

      WRITE(fp,*)
      WRITE(fp,'(A)') 'PLASMA DATA:'
      IF (nfluid.GT.0) THEN
        DO it = 1, ntube
          DO ion = 1, nion
            WRITE(tag,'(64X)')
            tag = tube_tag(tube(it)%type)
            IF (it.EQ.grid%isep ) tag = TRIM(tag)//', 1ST SEPARATRIX'
            IF (it.EQ.grid%isep2) tag = TRIM(tag)//', 2ND SEPARATRIX'

            WRITE(fp,*)
            WRITE(fp, 9) '   TUBE =',it,'    '//TRIM(tag)
            WRITE(fp,10) '   ION  =',ion
            WRITE(fp,'(3X,A8,3A11,3A10)') 
     .        'Cell','ne','ni','vi','M','Te','Ti'
            WRITE(fp,'(A8,3A11,3A10)') 
     .        '     ','(m-3)','(m-3)','(m s-1)',' ','(eV)','(eV)'
            WRITE(fp,'(11X,1P,3E11.2,0P,3F10.2)') 
     .        tube(it)%ne(LO),
     .        tube(it)%ni(LO,ion),
     .        tube(it)%vi(LO,ion),
     .       -tube(it)%machno(LO),
     .        tube(it)%te(LO),
     .        tube(it)%ti(LO,ion)
            i = 0
            DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)        
              i = i + 1
              cs = GetCs2(fluid(ic,ion)%te,fluid(ic,ion)%ti)
              WRITE(fp,'(I3,I8,1P,3E11.2,0P,3F10.2)') i,ic,
     .          fluid(ic,ion)%ne,
     .          fluid(ic,ion)%ni,
     .          fluid(ic,ion)%vi,
     .          fluid(ic,ion)%vi/cs,
     .          fluid(ic,ion)%te,
     .          fluid(ic,ion)%ti
            ENDDO
            WRITE(fp,'(11X,1P,3E11.2,0P,3F10.2)') 
     .        tube(it)%ne(HI),
     .        tube(it)%ni(HI,ion),
     .        tube(it)%vi(HI,ion),
     .        tube(it)%machno(HI),
     .        tube(it)%te(HI),
     .        tube(it)%ti(HI,ion)
          ENDDO
        ENDDO
      ENDIF

      WRITE(fp,*)
      WRITE(fp,'(A)') 'NEUTRAL DATA:'
      IF (npin.GT.0) THEN
        DO it = 1, ntube
          DO ion = 1, nion
            WRITE(fp,*)
            WRITE(fp,10) '   TUBE =',it
            WRITE(fp,10) '   ION  =',ion
            WRITE(fp,'(3X,A8,2(A11,6X),4A11)') 
     .        'Cell','n_atm','n_mol','ion','rec','Dalpha','Dgamma'
            WRITE(fp,'(11X,2(A11,6X),4A11)') 
     .        '(m-3)','(m-3)','(m-3 s-1)','(m-3 s-1)',
     .        '(m-3 s-1)','(m-3 s-1)'
            i = 0
            DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)        
              i = i +1
              WRITE(fp,'(I3,I8,2(1P,E11.2,0P,F6.2),1P,4E11.2,F6.2)') 
     .          i,ic,
     .          pin(ic,ion)%n_atm,
     .          pin(ic,ion)%n_atm / (fluid(ic,ion)%ne + EPS10),
     .          pin(ic,ion)%n_mol,
     .          pin(ic,ion)%n_mol / (fluid(ic,ion)%ne + EPS10),
     .          pin(ic,ion)%ion,
     .          pin(ic,ion)%rec,
     .          pin(ic,ion)%dalpha(1),
     .          pin(ic,ion)%dgamma(1),
     .          pin(ic,ion)%dgamma(1) / (pin(ic,ion)%dalpha(1) + EPS10)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      WRITE(fp,*)
      WRITE(fp,'(A)') 'FLUID GRID MAP DATA:'
      IF (nobj.GT.0) THEN
        DO it = 1, ntube
          WRITE(fp,*)
          WRITE(fp,10) '   TUBE =',it
          WRITE(fp,'(3X,6A8,2X,2A9)') 
     .      'Cell','Obj','Obj_12','Obj_23','Obj_34','Obj_41',
     .      'Tube_41','Tube_23'
          i = 0
          DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)        
            i = i + 1
            iobj = GetObject(ic,IND_CELL)
            WRITE(fp,'(I3,6I8,2X,2I9)') i,ic,iobj,
     .        obj(iobj)%omap(1:4),
     .        GetTube(obj(iobj)%omap(4),IND_OBJECT), 
     .        GetTube(obj(iobj)%omap(2),IND_OBJECT)
          ENDDO
        ENDDO      
      ENDIF


      WRITE(fp,*)
      WRITE(fp,'(A)') 'VECTOR FIELD DATA:'
      IF (nfield.GT.0) THEN
        DO it = 1, ntube
          WRITE(fp,*)
          WRITE(fp,10) '   TUBE =',it
          WRITE(fp,'(3X,A8,5A11)') 
     .      'Cell','Brat','Btot','Br','Bphi','Bz'
          WRITE(fp,'(11X,A8,5A11)') 
     .      '(Tesla)','(Tesla)','(Tesla)','(Tesla)','(Tesla)'
          i = 0
          DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)        
            i = i + 1
            WRITE(fp,'(I3,I8,5F11.5)') i,ic,
     .        field(ic)%bratio,
     .        field(ic)%b,
     .        field(ic)%br,
     .        field(ic)%bphi,
     .        field(ic)%bz
          ENDDO
        ENDDO
      ENDIF  


      WRITE(fp,*)
      WRITE(fp,'(A)') 'VERTEX DATA:'
      IF (nobj.GT.0.AND.nvtx.GT.0) THEN
        WRITE(fp,*)
        WRITE(fp,'(6X,A8,2A6,4I22)') 
     .    'Index','Tube','Cell',(i,i=1,4)
        WRITE(fp,'(6X,8X,12X,4A22)') 
     .    '(m)','(m)','(m)','(m)'
        DO iobj = 1, nobj
          IF (grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID.OR.  ! Only for fluid grid objects
     .        grp(obj(iobj)%group)%type  .NE.GRP_QUADRANGLE   ) CYCLE
          DO i = 1, obj(iobj)%nside
            CALL GetVertex(iobj,i,x(i),y(i))         
          ENDDO
          WRITE(fp,'(I6,I8,2I6,4(2X,2F10.6))') 
     .      iobj,
     .      -1,
     .      GetTube(iobj,IND_OBJECT),
     .      -1,
     .      (x(i),y(i),i=1,obj(iobj)%nside)
        ENDDO
      ENDIF  


      
      CLOSE(fp)
        
        
      RETURN
 98   CALL ER('DumpSolution','Trouble accessing file',*99)
 99   WRITE(0,*) '  FILE NAME= "'//TRIM(fname)//'"'
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE SaveFluidGridGeometry
      USE mod_interface
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_user
      IMPLICIT none


      INTEGER                 ion,igrp,iobj,iside,isrf,ivtx,it,ip
      INTEGER, ALLOCATABLE :: vlist(:)
      CHARACTER               dummy*1024


      WRITE(dummy,'(1024X)')

      ALLOCATE(vlist(nvtx))
      vlist = 0

      ion = 1

c...  ------------------------------------------------------------------
      CALL inOpenInterface('osm.idl.fluid_grid')

c      DO igrp = 1, ngrp
c        CALL inPutData(grp(igrp)%origin,'GRP_ORIGIN','none')
c        CALL inPutData(grp(igrp)%type  ,'GRP_TYPE'  ,'none')
c      ENDDO
 
      CALL inPutData(grid%isep,'GRD_ISEP' ,'none')
      CALL inPutData(grid%ipfz,'GRD_IPFZ' ,'none')

      CALL inPutData(ntube             ,'TUBE_N'   ,'none')
      CALL inPutData(tube(1:ntube)%psin,'TUBE_PSIN','none')
      CALL inPutData(tube(1:ntube)%rho ,'TUBE_RHO' ,'(m)')
      CALL inPutData(tube(1:ntube)%smax,'TUBE_L'   ,'(m)')

      DO iobj = 1, nobj
        igrp = obj(iobj)%group
        IF (grp(igrp)%origin.NE.GRP_MAGNETIC_GRID) CYCLE

        CALL inPutData(-1                       ,'OBJ_INDEX'   ,'none')
        CALL inPutData(obj(iobj)%index(IND_CELL),'OBJ_IND_CELL','none')
        DO it = 1, ntube   ! *** REPLACE WITH A FUNCTION ***
          IF (obj(iobj)%index(IND_CELL).GE.tube(it)%cell_index(LO).AND.
     .        obj(iobj)%index(IND_CELL).LE.tube(it)%cell_index(HI)) EXIT
        ENDDO
        IF (it.EQ.ntube+1)
     .    CALL ER('WriteFluidGridGeometry','TUBE index associated '//
     .            'with OBJECT not identified',*99)
        ip = obj(iobj)%index(IND_CELL) - tube(it)%cell_index(LO) + 1
        CALL inPutData(ip             ,'OBJ_IND_POS' ,'none')
        CALL inPutData(it             ,'OBJ_IND_TUBE','none')
        CALL inPutData(obj(iobj)%nside,'OBJ_NSIDE'   ,'none')

        DO iside = 1, obj(iobj)%nside
          WRITE(dummy,'(A,I1,500X)') 'OBJ_ISIDE_',iside
          CALL inPutData(obj(iobj)%iside(iside),TRIM(dummy),'none')
          WRITE(dummy,'(A,I1,500X)') 'OBJ_OMAP_' ,iside
          CALL inPutData(obj(iobj)%omap (iside),TRIM(dummy),'none')
        ENDDO
      ENDDO

      DO iobj = 1, nobj
        igrp = obj(iobj)%group
        IF (grp(igrp)%origin.NE.GRP_MAGNETIC_GRID) CYCLE
        DO iside = 1, obj(iobj)%nside
          isrf = obj(iobj)%iside(iside)
          IF (isrf.GT.0) THEN
            CALL inPutData(isrf             ,'SRF_INDEX' ,'none')
            CALL inPutData(srf(isrf)%type   ,'SRF_TYPE'  ,'none')            
            CALL inPutData(srf(isrf)%nvtx   ,'SRF_NVTX'  ,'none')            
            CALL inPutData(srf(isrf)%ivtx(1),'SRF_IVTX_1','none')  ! The only type of surface associated with 
            CALL inPutData(srf(isrf)%ivtx(2),'SRF_IVTX_2','none')  ! a fluid grid, at the moment...
            vlist(srf(isrf)%ivtx(1:2)) = srf(isrf)%ivtx(1:2)
          ENDIF
        ENDDO
      ENDDO

      DO ivtx = 1, nvtx
        IF (vlist(ivtx).EQ.0) CYCLE
        CALL inPutData(ivtx       ,'VTX_INDEX','none')
        CALL inPutData(vtx(1,ivtx),'VTX_1'    ,'m')
        CALL inPutData(vtx(2,ivtx),'VTX_2'    ,'m')
        CALL inPutData(vtx(3,ivtx),'VTX_3'    ,'m')
      ENDDO

      CALL inCloseInterface 

      DEALLOCATE(vlist)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE SaveWallGeometry
      USE mod_interface
      USE mod_sol28_wall
      IMPLICIT none

      INTEGER   iw
      CHARACTER dummy*1024

      WRITE(dummy,'(1024X)')

      CALL inOpenInterface('osm.idl.fluid_wall')
      DO iw = 1, nwall
        CALL inPutData(wall(iw)%class            ,'WALL_CLASS' ,'none')
        CALL inPutData(wall(iw)%index(WAL_GROUP ),'WALL_GROUP' ,'none')
        CALL inPutData(wall(iw)%index(WAL_INDEX ),'WALL_INDEX' ,'none')
        CALL inPutData(wall(iw)%index(WAL_TUBE  ),'WALL_TUBE'  ,'none')
        CALL inPutData(wall(iw)%index(WAL_TARGET),'WALL_TARGET','none')
        CALL inPutData(wall(iw)%v1(1),'WALL_V1_X','m')
        CALL inPutData(wall(iw)%v1(2),'WALL_V1_Y','m')
        CALL inPutData(wall(iw)%v2(1),'WALL_V2_X','m')
        CALL inPutData(wall(iw)%v2(2),'WALL_V2_Y','m')
      ENDDO
      CALL inCloseInterface 

      RETURN
 99   STOP
      END
c
c ======================================================================
c


