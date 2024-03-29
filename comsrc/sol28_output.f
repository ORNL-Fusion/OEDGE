c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE ListTargetData(fp,title)
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_sol28_targets
      IMPLICIT none      

      INTEGER  , INTENT(IN) :: fp
      CHARACTER, INTENT(IN) :: title*(*)

      REAL osm_GetGamma,osm_GetFlux,osm_GetHeatFlux

      INTEGER     itarget,itube,ion,i,ipos
      CHARACTER*2 target_tag(2) 
      REAL        area,heat_flux_density,heat(-1:ntarget)

      ion = 1

      target_tag(LO) = 'LO'
      target_tag(HI) = 'HI'

      IF (title.NE.'none') THEN 
        WRITE(fp,*)
        WRITE(fp,'(A)') 'TARGET DATA:'
        WRITE(fp,'(A)') '  '//TRIM(title)
      ENDIF

c     *** SOME INCONSISTENCY IN THIS ROUTINE WRT GAMMA, WHERE I'M PRESENTLY USING THE
c         (CRIPPLED) FORMULATION FROM STANGEBY, BUT IN SOME CASES I TRY TO BACK GAMMA 
c         OUT IN THE SOLVER, AND SO SHOULD BE USING THAT... ***

      DO itarget = 1, ntarget
        WRITE(fp,*)
        WRITE(fp,'(2X,A)') TRIM(target(itarget)%tag)
        WRITE(fp,'(A6,A8,A16,3A10,A6,4A10,2A9)') 
     .    'TUBE','psin','jsat','ne','ni','vi','M','pe','pi','Te','Ti',
     .    'gamma','heat'
        WRITE(fp,'(6X,8X,A16,3A10,6X,4A10,2A9)') 
     .    '(Amps)','(m-3)','(m-3)','(m s-1)','(?)','(?)',
     .    '(eV)','(eV)','(?)','(MW m-2)'
        DO i = 1, target(itarget)%nlist
          itube = target(itarget)%ilist(i)
          IF (tube(itube)%type.EQ.GRD_CORE) CYCLE
          ipos = target(itarget)%position
          area = 2.0 * V_PI * tube(itube)%rp (ipos) * 
     .                        tube(itube)%dds(ipos)
          heat_flux_density = osm_GetHeatFlux(ipos,itube) / area /1.0E+6
          WRITE(fp,'(I6,F8.4,1P,E16.6,3E10.2,0P,F6.2,
     .               1P,2E10.2,0P,2F10.4,F9.2,F9.2)')
     .      itube,tube(itube)%psin,
     .      tube(itube)%jsat  (ipos,ion),
     .      tube(itube)%ne    (ipos),
     .      tube(itube)%ni    (ipos,ion),
     .      tube(itube)%vi    (ipos,ion),
     .      tube(itube)%machno(ipos),
     .      tube(itube)%pe    (ipos),
     .      tube(itube)%pi    (ipos,ion),
     .      tube(itube)%te    (ipos),
     .      tube(itube)%ti    (ipos,ion),
     .      tube(itube)%gamma (ipos,ion),
     .      heat_flux_density
        ENDDO
      ENDDO

      heat = 0.0
      DO itarget = 1, ntarget
        ipos = target(itarget)%position
        DO i = 1, target(itarget)%nlist
          itube = target(itarget)%ilist(i)
          IF (tube(itube)%type.EQ.GRD_CORE) CYCLE
          heat(itarget) = heat(itarget) + osm_GetHeatFlux(ipos,itube)
        ENDDO
        heat(0) = heat(0) + heat(itarget)
      ENDDO
      heat = heat / 1.0E+6  ! Convert to MW m-2

      WRITE(fp,*)
      WRITE(fp,'(2X,A,F9.3,A)') 'TOTAL HEAT FLUX =',heat(0),' MW m-2'

      DO itarget = 1, ntarget
        WRITE(fp,*)
        WRITE(fp,'(2X,A,F7.3,A,F7.2,A)') 
     .    'HEAT FLUX TO '//TRIM(target(itarget)%tag)//
     .    ' TARGET =',heat(itarget),
     .    '  MW m-2',heat(itarget)/heat(0)*100.0, ' % '
        WRITE(fp,'(A6,2A8,2A12,2A10,A12)')
     .    'TUBE','psin','rho','gamma','isat','Te','Ti','heat flux'
        ipos = target(itarget)%position
        DO i = 1, target(itarget)%nlist
          itube = target(itarget)%ilist(i)
          IF (tube(itube)%type.EQ.GRD_CORE) CYCLE
          heat(-1) = osm_GetHeatFlux(ipos,itube) / 1.0E+6
          WRITE(fp,'(I6,2F8.4,1P,2E12.4,0P,2F10.2,
     .               1P,E12.4,0P,2(2X,F7.2,A))')
     .      itube,
     .      tube(itube)%psin,
     .      tube(itube)%rho ,
     .      osm_GetGamma   (ipos,itube),
     .      osm_GetFlux    (ipos,itube),
     .      tube(itube)%te (ipos),
     .      tube(itube)%ti (ipos,ion),
     .      heat(-1),heat(-1)/heat(itarget)*100.0,' %',
     .      heat(-1)/heat(0)*100.0,' %'
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
      SUBROUTINE GenerateOutputFiles(iitersol)
      USE mod_interface
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_sol28_targets
      USE mod_user
      IMPLICIT none
  
      INTEGER, INTENT(IN) :: iitersol

      REAL GetCs2,GetCellPressure

      INTEGER     ion,itube,ipos,icell,itar,ic1,ic2,itube1,itube2,
     .            location,i1,ival(10),itarget
      CHARACTER*7 tag
      CHARACTER*2 target_tag(2)
      CHARACTER*4 iter
      REAL        val(10),rdum

      INTEGER it(ncell),ic(ncell)
      REAL    cs(ncell),M(ncell),pe(ncell),pi(ncell)  

      LOGICAL first_call
      DATA    first_call /.TRUE./
      SAVE


      CALL DumpData_OSM('output.end','Simulation complete')

c...  Save solution:
      CALL SaveGrid('osm.raw')

      CALL SaveFluidGridGeometry
      CALL SaveWallGeometry

      IF (iitersol.EQ.-999) THEN 
        iter = '    '
      ELSE
        WRITE(iter,'(A,I3.3)') '_',iitersol 
      ENDIF

c...  Output files for CORTEX:

      ion = 1

      itube1 = 1
      itube2 = ntube
c      itube1 = 4
c      itube2 = 4  ! ntube

c...  ------------------------------------------------------------------
      CALL inOpenInterface('osm.idl.fluid_plasma'//iter,ITF_WRITE)

      DO itube = itube1, itube2
        ic1 = tube(itube)%cell_index(LO)
        ic2 = tube(itube)%cell_index(HI)
        DO icell = ic1, ic2
          it (icell) = itube
          ic (icell) = icell
          cs (icell) = GetCs2(fluid(icell,ion)%te,fluid(icell,ion)%ti)
          M  (icell) = fluid(icell,ion)%vi / cs(icell)
          pe (icell) = GetCellPressure(icell,1)
          pi (icell) = GetCellPressure(icell,2)
        ENDDO

        CALL inPutData(tube(itube)%psin,'PSIN','none')
        CALL inPutData(tube(itube)%rho ,'RHO' ,'none')

        CALL inPutData(it(ic1:ic2),'TUBE' ,'none')
        CALL inPutData(ic(ic1:ic2),'INDEX','none')
        CALL inPutData(cell(ic1:ic2)%s,'S','m')
        CALL inPutData(cell(ic1:ic2)%p,'P','m')
        CALL inPutData(cell(ic1:ic2)%cencar(1),'R','m')
        CALL inPutData(cell(ic1:ic2)%cencar(2),'Z','m')
        CALL inPutData(cell(ic1:ic2)%sbnd(1),'SBND1','m')
        CALL inPutData(cell(ic1:ic2)%sbnd(2),'SBND2','m')
        CALL inPutData(cell(ic1:ic2)%pbnd(1),'PBND1','m')
        CALL inPutData(cell(ic1:ic2)%pbnd(2),'PBND2','m')
        CALL inPutData(fluid(ic1:ic2,ion)%ne,'NE'    ,'m-3')        
        CALL inPutData(fluid(ic1:ic2,ion)%vi,'VI'    ,'m s-1')        
        CALL inPutData(cs    (ic1:ic2)      ,'CS'    ,'m s-1')
        CALL inPutData(M     (ic1:ic2)      ,'MACHNO','m s-1')
        CALL inPutData(pe    (ic1:ic2)      ,'PE'    ,'Pa')
        CALL inPutData(pi    (ic1:ic2)      ,'PI'    ,'Pa')
        CALL inPutData(fluid(ic1:ic2,ion)%te,'TE'    ,'eV')
        CALL inPutData(fluid(ic1:ic2,ion)%ti,'TI'    ,'eV')
      ENDDO
      CALL inCloseInterface 

c...  ------------------------------------------------------------------
      CALL inOpenInterface('osm.idl.fluid_sources'//iter,ITF_WRITE)
      DO itube = itube1, itube2
        ic1 = tube(itube)%cell_index(LO)
        ic2 = tube(itube)%cell_index(HI)
c....   Cells centres:
        CALL inPutData(it(ic1:ic2),'TUBE' ,'none')
        CALL inPutData(ic(ic1:ic2),'INDEX','none')
        CALL inPutData(cell(ic1:ic2)%s,'S','m')
        CALL inPutData(cell(ic1:ic2)%p,'P','m')
        CALL inPutData(cell(ic1:ic2)%cencar(1),'R','m')
        CALL inPutData(cell(ic1:ic2)%cencar(2),'Z','m')
        CALL inPutData(fluid(ic1:ic2,ion)%parsrc,'PAR_NET','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%parion,'PAR_ION','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%parrec,'PAR_REC','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%parano,'PAR_ANO','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%parusr,'PAR_USR','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%momsrc,'MOM_NET','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%momvol,'MOM_VOL','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%momano,'MOM_ANO','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%momusr,'MOM_USR','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%enesrc,'ENE_NET','?')
        CALL inPutData(fluid(ic1:ic2,ion)%eneion,'ENE_ION','?')
        CALL inPutData(fluid(ic1:ic2,ion)%enerec,'ENE_REC','?')
        CALL inPutData(fluid(ic1:ic2,ion)%eneano,'ENE_FIT','?')
        CALL inPutData(fluid(ic1:ic2,ion)%eneusr,'ENE_USR','?')
        CALL inPutData(fluid(ic1:ic2,ion)%enisrc,'ENI_NET','?')
        CALL inPutData(fluid(ic1:ic2,ion)%eniion,'ENI_ION','?')
        CALL inPutData(fluid(ic1:ic2,ion)%eniano,'ENI_FIT','?')
        CALL inPutData(fluid(ic1:ic2,ion)%eniusr,'ENI_USR','?')
      ENDDO
      CALL inCloseInterface 

c...  ------------------------------------------------------------------
      IF (ALLOCATED(pin)) THEN
        CALL inOpenInterface('osm.idl.fluid_eirene'//iter,ITF_WRITE)
        DO itube = itube1, itube2
          ic1 = tube(itube)%cell_index(LO)
          ic2 = tube(itube)%cell_index(HI)
c....     Cells centres:
          CALL inPutData(ic(ic1:ic2),'INDEX','none')
          CALL inPutData(ic(ic1:ic2),'POS','none')
          CALL inPutData(it(ic1:ic2),'TUBE' ,'none')
          CALL inPutData(cell(ic1:ic2)%s,'S','m')
          CALL inPutData(cell(ic1:ic2)%cencar(1),'R','m')
          CALL inPutData(cell(ic1:ic2)%cencar(2),'Z','m')
          CALL inPutData(pin(ic1:ic2,ion)%ion  ,'ION_NET','?')        
          CALL inPutData(pin(ic1:ic2,ion)%rec  ,'REC_NET','?')        
          CALL inPutData(pin(ic1:ic2,ion)%mom  ,'MOM_NET','?')
          CALL inPutData(pin(ic1:ic2,ion)%qe   ,'QE_NET' ,'?')
          CALL inPutData(pin(ic1:ic2,ion)%qi   ,'QI_NET' ,'?')
          CALL inPutData(pin(ic1:ic2,ion)%n_atm,'ATM_DENS','?')
          CALL inPutData(pin(ic1:ic2,ion)%n_mol,'MOL_DENS','?')
          CALL inPutData(pin(ic1:ic2,ion)%dalpha(1),'BALMER_ALPHA','?')
          CALL inPutData(pin(ic1:ic2,ion)%dgamma(1),'BALMER_GAMMA','?')
        ENDDO
        CALL inCloseInterface 
      ENDIF
c...  ------------------------------------------------------------------
      IF (ALLOCATED(field)) THEN
        CALL inOpenInterface('osm.idl.fluid_fields'//iter,ITF_WRITE)
        DO itube = itube1, itube2
          ic1 = tube(itube)%cell_index(LO)
          ic2 = tube(itube)%cell_index(HI)
c....     Cells centres:
          CALL inPutData(ic(ic1:ic2),'INDEX','none')
          CALL inPutData(ic(ic1:ic2),'POS','none')
          CALL inPutData(it(ic1:ic2),'TUBE' ,'none')
          CALL inPutData(cell(ic1:ic2)%s,'S','m')
          CALL inPutData(field(ic1:ic2)%epot  ,'EPOT'  ,'?')        
          CALL inPutData(field(ic1:ic2)%efield,'EFIELD','?')        
        ENDDO
        CALL inCloseInterface 
      ENDIF
c...  ------------------------------------------------------------------
      CALL inOpenInterface('osm.idl.params'//iter,ITF_WRITE)
      CALL inPutData(2.0,'flupar mass','amu')
      CALL inCloseInterface 
c
c     ------------------------------------------------------------------
c     Write out target data:
c
c     If changing anything here, need to change it in DumpDataToIDL in
c     SLoutplot.f as well, so that the OSM and OUT generated 
c     idl.fluid_targets files remain in sync.
c
      IF (.NOT.ALLOCATED(target)) GOTO 20

      CALL inOpenInterface('osm.idl.fluid_targets'//iter,ITF_WRITE)

      target_tag(LO) = 'LO'
      target_tag(HI) = 'HI'

      DO itube = 1, ntube
        IF (itube.LT.grid%isep) CYCLE
        CALL inPutData(itube           ,'TAR_TUBE','none')                    
        CALL inPutData(tube(itube)%ir  ,'TAR_RING','none')                    
        CALL inPutData(tube(itube)%psin,'TAR_PSIN','none')                    
        CALL inPutData(tube(itube)%rho ,'TAR_RHO' ,'m') 
        DO ipos = LO, HI  
         DO itar = 1, ntarget
           DO i1 = 1, target(itar)%nlist
             IF (itube.EQ.target(itar)%ilist(i1).AND.
     .           ipos .EQ.target(itar)%position) EXIT
           ENDDO
           IF (i1.NE.target(itar)%nlist+1) EXIT
         ENDDO
         IF (itar.EQ.ntarget+1.AND.i1.EQ.target(ntarget)%nlist+1) THEN
           IF (first_call) THEN
             CALL WN('GenerateOutputFiles','Some target group '//
     .               'identifiers not found')
             first_call = .FALSE.
           ENDIF
           itarget  = -1
           location = -1
         ELSE
           itarget  = itar
           location = target(itar)%location
         ENDIF
         tag = 'TAR_'//target_tag(ipos)//'_'
         IF (ipos.EQ.LO) THEN
           CALL inPutData(0.0             ,tag//'S','m')                    
           CALL inPutData(0.0             ,tag//'P','m') ! Should have THETA here as well, and the target extent and orientation to the indicent field line               
         ELSE
           CALL inPutData(tube(itube)%smax,tag//'S','m')                    
           CALL inPutData(tube(itube)%pmax,tag//'P','m')                    
         ENDIF
c         CALL inPutData(tube(itube)%rp(ipos),tag//'R','m')                    
c         CALL inPutData(tube(itube)%zp(ipos),tag//'Z','m') 
         CALL inPutData(itarget ,tag//'TARGET_INDEX','none')                    
         CALL inPutData(location,tag//'LOCATION'    ,'none')                    
         CALL inPutData(tube(itube)%jsat(ipos,ion),tag//'JSAT','Amps')                    
         CALL inPutData(tube(itube)%ne  (ipos    ),tag//'NE'  ,'m-3')                    
         CALL inPutData(tube(itube)%vi  (ipos,ion),tag//'V'   ,'m s-1')                    
         CALL inPutData(tube(itube)%machno(ipos)  ,tag//'MACHNO','none')                    
         CALL inPutData(tube(itube)%pe  (ipos    ),tag//'PE'  ,'m-3 eV')                    
         CALL inPutData(tube(itube)%pi  (ipos,ion),tag//'PI'  ,'m-3 eV')                    
         CALL inPutData(tube(itube)%te  (ipos    ),tag//'TE'  ,'eV')                    
         CALL inPutData(tube(itube)%ti  (ipos,ion),tag//'TI'  ,'eV')                    
        ENDDO
      ENDDO
      CALL inCloseInterface

 20   CONTINUE

c...  ------------------------------------------------------------------
      CALL inOpenInterface('osm.idl.osm_nodes'//iter,ITF_WRITE)
      DO itube = 1, ntube
        CALL inPutData(itube             ,'TUBE'  ,'N/A')
        CALL inPutData(store_sopt (itube),'S_OPT' ,'N/A')
        CALL inPutData(store_mnode(itube),'M_NODE','N/A')
        CALL inPutData(store_nnode(itube),'N_NODE','N/A')
        DO i1 = 1, MAXVAL(store_nnode(1:ntube))
          WRITE(tag,'(A,I1,A)') 'NODE_',i1,'_'
          CALL inPutData(store_node(i1,itube)%s,tag//'S','m')        
          CALL inPutData(store_node(i1,itube)%jsat(ion),tag//'JSAT','A')        
          CALL inPutData(store_node(i1,itube)%ne,tag//'DENS','m-3')        
          CALL inPutData(store_node(i1,itube)%pe,tag//'PE','?')        
          CALL inPutData(store_node(i1,itube)%te,tag//'TE','eV')        
          CALL inPutData(store_node(i1,itube)%ti(ion),tag//'TI','eV')        
        ENDDO
      ENDDO
      DO i1 = 1, osmnnode
        CALL inPutData(osmnode(i1)%type         ,'TYPE'       ,'N/A')
        CALL inPutData(osmnode(i1)%tube_range(1),'TUBE_RANGE1','N/A')
        CALL inPutData(osmnode(i1)%tube_range(2),'TUBE_RANGE2','N/A')
        CALL inPutData(osmnode(i1)%rad_x        ,'RAD_X','(m)')
        CALL inPutData(osmnode(i1)%rad_y        ,'RAD_Y','(m)')
      ENDDO
      CALL inCloseInterface 


c...  Give the user a chance to dump some output:
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
c            WRITE(fp,*) 'GETOBJECT: ',i1,iobj,imap,idim2
            IF (imap.EQ.0) CYCLE
            IF (i1.EQ.IND_IK.OR.i1.EQ.IND_IR.OR.i1.EQ.IND_IS) CYCLE
c            WRITE(fp,*) '         : ',iobj,i1,imap,obj_index_map(imap,i1)
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
      REAL*8    x(10),y(10),zp(2),x1,x2,y1,y2

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
        WRITE(fp,'(2A6,A5,2A9,11A10)') 
     .    'Index','Type','#','Cell_LO','Cell_HI',
     .    'psi_n','rho(m)','L (m)',
     .    'costet_LO','costet_HI',
     .    'rp_LO'    ,'rp_HI'    ,
     .    'zp_LO'    ,'zp_HI'    ,
     .    'dds_LO'   ,'dds_HI'   
        DO it = 1, ntube

          iobj = GetObject(tube(it)%cell_index(LO),IND_CELL)
          CALL GetVertex(iobj,1,x1,y1)         
          CALL GetVertex(iobj,2,x2,y2)         
          zp(LO) = 0.5D0 * (y1 + y2)
          iobj = GetObject(tube(it)%cell_index(HI),IND_CELL)
          CALL GetVertex(iobj,3,x1,y1)         
          CALL GetVertex(iobj,4,x2,y2)         
          zp(HI) = 0.5D0 * (y1 + y2)

          WRITE(fp,'(I6,A6,I5,2I9,13F10.4)')
     .      it,
     .      tube_tag(tube(it)%type),
     .      tube(it)%cell_index(HI)-tube(it)%cell_index(LO)+1,
c     .      tube(it)%n,
     .      tube(it)%cell_index(LO:HI),
     .      tube(it)%psin,
     .      tube(it)%rho,
     .      tube(it)%smax,
     .      tube(it)%costet(1:2),
     .      tube(it)%rp    (1:2),
     .      zp(1:2),
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

      CALL ListTargetData(fp,'none')

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
            WRITE(fp,'(3X,A8,10A11)') 
     .        'Cell','bratio','s','sbnd1','sbnd2','p','q','R','Z','vol',
     .        'metric'
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
              WRITE(fp,'(I3,I8,6F11.5,4F11.5)') i,ic,
     .          field(ic)%bratio,
     .          cell (ic)%s,
     .          cell (ic)%sbnd(1:2),
     .          cell (ic)%p,
     .          -1.0,
     .          cell (ic)%cencar(1:2),
     .          cell (ic)%vol,
     .          cell (ic)%metric
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
     .          fluid(ic,ion)%vi/(cs+eps10),
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
      CALL inOpenInterface('osm.idl.fluid_grid',ITF_WRITE)

c      DO igrp = 1, ngrp
c        CALL inPutData(grp(igrp)%origin,'GRP_ORIGIN','none')
c        CALL inPutData(grp(igrp)%type  ,'GRP_TYPE'  ,'none')
c      ENDDO
c      write(0,*) 'isep' ,grid%isep 
      CALL inPutData(grid%isep,'GRD_ISEP' ,'none')
c      write(0,*) 'ipfz' ,grid%ipfz
      CALL inPutData(grid%ipfz,'GRD_IPFZ' ,'none')

      CALL inPutData(grid%r0,'GRD_R0','m')
      CALL inPutData(grid%z0,'GRD_Z0','m')

      CALL inPutData(grid%rxpt,'GRD_RXPT','m')
      CALL inPutData(grid%zxpt,'GRD_ZXPT','m')

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

      IF (nwall.EQ.0) RETURN

      CALL inOpenInterface('osm.idl.fluid_wall',ITF_WRITE)
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


