c     -*-Fortran-*-
c
c ======================================================================
c
c input.f
c
c ======================================================================
c
      SUBROUTINE User_InitializeOptions
      USE mod_sol28_global
      USE mod_user
      IMPLICIT none

      opt_user%u_mom = 0

      


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_AllocateArrays
      USE mod_sol28_global
      USE mod_user
      IMPLICIT none

c      WRITE(0,*) 'NCELL IN SETUPARRAYS:',ncell
      ALLOCATE(stuff(ncell))      


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_DeallocateArrays
      USE mod_sol28_global
      USE mod_user
      IMPLICIT none


      DEALLOCATE(stuff)      


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_LoadOptions(fp,itag,buffer)
      USE mod_sol28_global
      USE mod_user
      IMPLICIT none

      INTEGER   fp,itag
      CHARACTER buffer*(*)

      SELECTCASE (buffer(3:itag-1))

        CASE('S U_MOM')
          CALL ReadOptionI(buffer,2,opt_user%u_mom)


        CASE DEFAULT
          CALL ER('User_LoadOptions','Unrecognized input tag',*99)
      ENDSELECT 

      RETURN
 99   WRITE(0,*) 'TAG: >'//buffer(3:itag-1)//'<'
      STOP
      END
c
c ======================================================================
c
c output.f
c
c ======================================================================
c
      SUBROUTINE User_GenerateOutputFiles
      USE mod_interface
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_sol28_targets
      USE mod_user
      IMPLICIT none

      REAL GetCs2,GetCellPressure

      INTEGER     ion,itube,ipos,icell,itar,ic1,ic2,itube1,itube2,
     .            location,i1,ival(10),itarget
      CHARACTER*7 tag
      CHARACTER*2 target_tag(2)
      REAL        val(10),rdum

      INTEGER it(ncell),ic(ncell)
      REAL    cs(ncell),M(ncell),pe(ncell),pi(ncell)  

      LOGICAL first_call
      DATA    first_call /.TRUE./
      SAVE

      ion = 1

      itube1 = 1
      itube2 = ntube
c      itube1 = 4
c      itube2 = 4  ! ntube

c...  ------------------------------------------------------------------
      CALL inOpenInterface('osm.idl.fluid_plasma',ITF_WRITE)

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
        CALL inPutData(cell(ic1:ic2)%sbnd(1),'SBND1' ,'m')
        CALL inPutData(cell(ic1:ic2)%sbnd(2),'SBND2' ,'m')
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
      CALL inOpenInterface('osm.idl.fluid_sources',ITF_WRITE)
      DO itube = itube1, itube2
        ic1 = tube(itube)%cell_index(LO)
        ic2 = tube(itube)%cell_index(HI)
c....   Cells centres:
        CALL inPutData(it(ic1:ic2),'TUBE' ,'none')
        CALL inPutData(ic(ic1:ic2),'INDEX','none')
        CALL inPutData(cell(ic1:ic2)%s,'S','m')
        CALL inPutData(fluid(ic1:ic2,ion)%parsrc,'PAR_NET','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%parion,'PAR_ION','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%parrec,'PAR_REC','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%parano,'PAR_ANO','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%parusr,'PAR_USR','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%momsrc,'MOM_NET','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%momvol,'MOM_VOL','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%momano,'MOM_ANO','?')        
        CALL inPutData(fluid(ic1:ic2,ion)%enesrc,'ENE_NET','?')
        CALL inPutData(fluid(ic1:ic2,ion)%eneion,'ENE_ION','?')
        CALL inPutData(fluid(ic1:ic2,ion)%enerec,'ENE_REC','?')
        CALL inPutData(fluid(ic1:ic2,ion)%eneano,'ENE_FIT','?')
        CALL inPutData(fluid(ic1:ic2,ion)%eneusr,'ENE_USR','?')
        CALL inPutData(fluid(ic1:ic2,ion)%enisrc,'ENI_NET','?')
      ENDDO
      CALL inCloseInterface 

c...  ------------------------------------------------------------------
      IF (ALLOCATED(pin)) THEN
        CALL inOpenInterface('osm.idl.fluid_eirene',ITF_WRITE)
        DO itube = itube1, itube2
          ic1 = tube(itube)%cell_index(LO)
          ic2 = tube(itube)%cell_index(HI)
c....     Cells centres:
          CALL inPutData(ic(ic1:ic2),'INDEX','none')
          CALL inPutData(ic(ic1:ic2),'POS','none')
          CALL inPutData(it(ic1:ic2),'TUBE' ,'none')
          CALL inPutData(cell(ic1:ic2)%s,'S','m')
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
      CALL inOpenInterface('osm.idl.params',ITF_WRITE)
      CALL inPutData(2.0,'flupar mass','amu')
      CALL inCloseInterface 
c
c     ------------------------------------------------------------------
c     Write out target data:
c
      IF (.NOT.ALLOCATED(target)) GOTO 20

      CALL inOpenInterface('osm.idl.fluid_targets',ITF_WRITE)
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
             CALL WN('User_GenerateOutputFiles','Some target group '//
     .               'identifiers not found')
             first_call = .FALSE.
           ENDIF
           CYCLE
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
      CALL inOpenInterface('osm.idl.osm_nodes',ITF_WRITE)
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

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c source.f
c
c ======================================================================
c
      SUBROUTINE User_VolumeParIonSource(target,source)
      USE mod_sol28_solver
      USE mod_user
      IMPLICIT none
 
      INTEGER,INTENT(IN)  :: target
      REAL*8 ,INTENT(OUT) :: source(*)

      INTEGER ic1,ic2

      ic1 = icbnd1(target)
      ic2 = icbnd2(target)

      source(ic1:ic2) = 0.0D0

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_VolumeParRecSource(target,source)
      USE mod_sol28_solver
      USE mod_user
      IMPLICIT none
 
      INTEGER,INTENT(IN)  :: target
      REAL*8 ,INTENT(OUT) :: source(*)

      INTEGER ic1,ic2

      ic1 = icbnd1(target)
      ic2 = icbnd2(target)

      source(ic1:ic2) = 0.0D0

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_VolumeParAnoSource(target,source)
      USE mod_sol28_solver
      USE mod_user
      IMPLICIT none
 
      INTEGER,INTENT(IN)  :: target
      REAL*8 ,INTENT(OUT) :: source(*)

      INTEGER ic1,ic2

      ic1 = icbnd1(target)
      ic2 = icbnd2(target)

      source(ic1:ic2) = 0.0D0

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_ParticleSource(ion,target,source)
      USE mod_sol28_solver
      USE mod_user
      IMPLICIT none
 
      INTEGER,INTENT(IN)  :: ion,target
      REAL*8 ,INTENT(OUT) :: source(*)

      INTEGER ic1,ic2

      ic1 = icbnd1(target)
      ic2 = icbnd2(target)

      source(ic1:ic2) = 0.0D0

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_MomentumSource(ion,target,source)
      USE mod_sol28_solver
      USE mod_user
      IMPLICIT none
 
      INTEGER,INTENT(IN)  :: ion,target
      REAL*8 ,INTENT(OUT) :: source(*)

      INTEGER ic1,ic2

      ic1 = icbnd1(target)  ! These set the cell index bounds for half flux-tube...
      ic2 = icbnd2(target)

c      STOP 'HOLY CRAP'

      ! DO ic = ic1, ic2
      !   R = cell(ic)%cencar(1)
      !   Z = cell(ic)%cencar(2)
      !   S = cell(ic)%s
      !   V = cell(ic)%vol
      ! ENDDO  

      source(ic1:ic2) = 0.0D0

      !ico = tube%cell_index(1) - 1
      ! stuff(ic1+ico:ic2+ico) = source(ic1:ic2)

      RETURN
 99   STOP
      END
c ======================================================================
c
      SUBROUTINE User_VolumeEneRecSource(target,source)
      USE mod_sol28_solver
      USE mod_user
      IMPLICIT none
 
      INTEGER,INTENT(IN)  :: target
      REAL*8 ,INTENT(OUT) :: source(*)

      INTEGER ic1,ic2

      ic1 = icbnd1(target)
      ic2 = icbnd2(target)

      source(ic1:ic2) = 0.0D0

      RETURN
 99   STOP
      END
c ======================================================================
c
      SUBROUTINE User_VolumeEneIonSource(target,source)
      USE mod_sol28_solver
      USE mod_user
      IMPLICIT none
 
      INTEGER,INTENT(IN)  :: target
      REAL*8 ,INTENT(OUT) :: source(*)

      INTEGER ic1,ic2

      ic1 = icbnd1(target)
      ic2 = icbnd2(target)

      source(ic1:ic2) = 0.0D0

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_DumpData(fp,fname,title)
      USE mod_sol28_global
      USE mod_user
      IMPLICIT none

      INTEGER       fp
      CHARACTER*(*) fname,title 

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c main.f
c
c ======================================================================
c
      SUBROUTINE User_SetupSolverOptions(count)
      USE mod_sol28_solver

      INTEGER, INTENT(IN):: count

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE User_MainLoop(condition,count)
      USE mod_sol28_global
      USE mod_user
      IMPLICIT none

      LOGICAL             :: condition
      INTEGER, INTENT(IN) :: count
      
      RETURN
 99   STOP
      END
