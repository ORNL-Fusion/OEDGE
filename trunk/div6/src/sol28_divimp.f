c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE SetupSOL28
      USE mod_sol28
      USE mod_sol28_global
      USE mod_legacy
      USE mod_geometry
      IMPLICIT none

c      RETURN

      CALL MapRingstoTubes

      CALL SaveGeometryData('osm_geometry.raw')

      CALL LoadLegacyData('osm_legacy.raw')

      CALL SetTargetConditions

      IF (opt%osm_load.NE.0) CALL LoadReferenceSolution(1)

c...  Clear geometry arrays:
      nobj = 0
      nsrf = 0
      nvtx = 0
      IF (ALLOCATED(obj)) DEALLOCATE(obj) 
      IF (ALLOCATED(srf)) DEALLOCATE(srf) 
      IF (ALLOCATED(vtx)) DEALLOCATE(vtx)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE ExecuteSOL28(irstart,irend,ikopt,sloutput)
      USE mod_sol28
      USE mod_sol28_global
      USE mod_geometry
      USE mod_legacy
      IMPLICIT none

      INTEGER irstart,irend,ikopt
      LOGICAL sloutput

      INTEGER itube,itube1,itube2,status
      LOGICAL cont

c...  Need to reload geometry data in case Eirene has been called (hopefully
c     this won't be required in the future...):
      CALL LoadObjects('osm_geometry.raw',status)
      IF (status.NE.0) CALL ER('ExecuteSOL28','Unable to load '//
     .                         'geometry data',*99)

c...  Load up PIN data if available:
      IF (opt%pin_data) CALL MapNeutralstoTubes

c      CALL MapRingstoTubes
c      CALL LoadLegacyData('osm_legacy.raw')
c      CALL SetTargetConditions

      itube1 = 0
      itube2 = 0
      DO itube = 1, ntube
        IF (tube(itube)%ir.EQ.irstart) itube1 = itube
        IF (tube(itube)%ir.EQ.irend  ) itube2 = itube
      ENDDO
      IF (itube2.EQ.0.AND.irend.EQ.ntube+2) itube2 = ntube

c...  Call SOL28 plasma solver:
      CALL MainLoop(itube1,itube2,ikopt,sloutput)

c...  Fill DIVIMP arrays:
      CALL MapTubestoRings(irstart,irend)

c...  Clear geometry arrays:
      nobj = 0
      nsrf = 0
      nvtx = 0
      IF (ALLOCATED(obj)) DEALLOCATE(obj) 
      IF (ALLOCATED(srf)) DEALLOCATE(srf) 
      IF (ALLOCATED(vtx)) DEALLOCATE(vtx)

c      CALL CleanUp

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE CloseSOL28
      IMPLICIT none

c      RETURN

c...  Save solution:
      CALL SaveGrid('osm.raw')

c...  Clear memory:
      CALL CleanUp

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE AssignNodeValues(itube,nnode,mnode,node)
      USE mod_sol28
      USE mod_sol28_global
      IMPLICIT none

      INTEGER itube,nnode,mnode       
      TYPE(type_node) :: node(*)

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      REAL GetJsat

      INTEGER ik,ir,i1,id,id1,id7
      LOGICAL new
      REAL    te(0:6),ne(0:6),nf,s(0:6)

      STOP 'ROUTINE TAGGED FOR DELETION'
    

      ir = tube(itube)%ir

      new = .TRUE.

c     get rid of do loop, put job of finding where in the SOL28 input listing should be used into a subroutine

      CALL FindS28Parameters_V3(ir,te,ne,nf,s,new)



      IF (.TRUE.) THEN

        s(0) = 0.0
        s(6) = ksmaxs(ir)

        DO i1 = 0, 6
          DO ik = 1, nks(ir)
            IF (s(i1).GE.ksb(ik-1,ir).AND.s(i1).LE.ksb(ik,ir)) THEN
              node(i1+1)%icell = ik
              IF (i1.EQ.3) s(3) = kss(ik,ir)
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

        id1 = idds(ir,2)
        id7 = idds(ir,1)
        IF (node(1)%ne.EQ.0.0) node(1)%ne = knds (id1)
        IF (node(7)%ne.EQ.0.0) node(7)%ne = knds (id7)
        IF (node(1)%te.EQ.0.0) node(1)%te = kteds(id1)
        IF (node(7)%te.EQ.0.0) node(7)%te = kteds(id7)
        IF (node(1)%ti(1).EQ.0.0) node(1)%ti(1) = ktids(id1)
        IF (node(7)%ti(1).EQ.0.0) node(7)%ti(1) = ktids(id7)

c REAL FUNCTION GetJsat(te,ti,ne,v)
        id = id1
        node(1)%jsat(1)= GetJsat(kteds(id),ktids(id),knds(id),kvds(id))
        node(1)%jsat(1) = -ABS(node(1)%jsat(1))
        id = id7
        node(7)%jsat(1)= GetJsat(kteds(id),ktids(id),knds(id),kvds(id))
        node(7)%jsat(1) = -ABS(node(7)%jsat(1))
      ENDIF


      RETURN
99    CONTINUE
      WRITE(0,*) 'IK,IR=',ik,ir,i1,osms28(i1,1)
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE MapNeutralstoTubes
      USE mod_sol28_global
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ion,ncell1,ike,ir,cind1,cind2

c...  Copy PIN data:
      ion = 1
      ncell1 = 0
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1

        cind1 = ncell1 + 1
        cind2 = ncell1 + ike 

        pin(cind1:cind2,ion)%ion = pinion(1:ike,ir)
        pin(cind1:cind2,ion)%rec = pinrec(1:ike,ir)
        pin(cind1:cind2,ion)%mom = pinmp (1:ike,ir)

        ncell1 = cind2
      ENDDO

      IF (ncell.NE.ncell1) 
     .  CALL ER('MapNeutraltoTubes','NCELL1.NE.NCELL',*99)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE MapRingstoTubes
      USE mod_sol28_global
      USE mod_geometry
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ik,ike,ir,id,cind1,cind2,i1,i2,iside,ik1,ir1,ir2,ion,iobj,
     .        fp,idum1
      LOGICAL pfz_ring,load_bfield_data
      REAL    rdum1
      REAL, ALLOCATABLE:: bfield_data(:,:,:)
      REAL*8  a(3)
      TYPE(type_srf   ) newsrf
      TYPE(type_object) newobj


      load_bfield_data = .FALSE.

      IF (load_bfield_data) THEN
        ALLOCATE(bfield_data(MAXNKS,MAXNRS,5))
        fp = 99
        OPEN(UNIT=fp,FILE='objects.bfield',ACCESS='SEQUENTIAL',
     .       STATUS='OLD',ERR=95)      
        READ(fp,*)
        READ(fp,*)
        DO ir = 1, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          DO ik = 1, nks(ir)
            READ(fp,*) (rdum1,i1=1,4),ik1,ir1,rdum1,
     .                 bfield_data(ik,ir,1:5)
            IF (ik1.NE.ik.OR.ir1.NE.ir) 
     .        CALL ER('MapRingsToTubes','Index mismatch',*99)
          ENDDO
        ENDDO
        CLOSE(fp)
      ENDIF


c *CRUDE* Should be elsewhere... can IDRING take the place of RINGTYPE for EIRENE04? 
      ringtype = 0
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) THEN
          ringtype(ir) = BOUNDARY
          CYCLE
        ENDIF
        IF (ir.LT.irsep) ringtype(ir) = CORE
        IF (ir.GE.irsep.AND.ir.LT.irwall) ringtype(ir) = SOL1
        IF (ir.GT.irtrap.AND.ir.LE.nrs) ringtype(ir) = PFZ
      ENDDO
c...  Special for secondary PFZ in generalized grids:
      IF (grdnmod.NE.0) THEN
        DO ik = 1, nks(irwall)
c...      Scan over IRWALL and check all rings that point to IRWALL to
c         see if *all* of that ring is pointing at IRWALL, in which
c         case it must be a PFZ ring (could just be a discontinuous
c         SOL target otherwise):
          IF (irins(ikins(ik,irwall),irins(ik,irwall)).EQ.irwall.AND.
     .        ringtype(irins(ik,irwall)).NE.PFZ) THEN
            ir1 = irins(ik,irwall)  
            pfz_ring = .TRUE.
            DO ik1 = 1, nks(ir1)
              IF (irins(ik1,ir1).NE.irwall) pfz_ring = .FALSE.
            ENDDO
            IF (pfz_ring) THEN
              DO ir2 = ir1, nrs
                ringtype(ir2) = PFZ
                IF     (irouts(1,ir2).EQ.irouts(nks(ir2),ir2))THEN
                ELSEIF (irouts(1       ,ir2).EQ.irwall.OR.
     .                  irouts(nks(ir2),ir2).EQ.irwall) THEN
c...              Discontinuous target!
                  STOP 'THIS IS A PROBLEM'
                ELSE
                  EXIT
                ENDIF
              ENDDO
            ENDIF
          ENDIF

        ENDDO
      ENDIF

c...  Count number of cells:
      ntube = 0
      ncell = 0
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ntube = ntube + 1
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1
        ncell = ncell + ike 
      ENDDO

c...  Declare global arrays:
      nfield    = ncell
      npin      = ncell
      nphoton   = 1
      ndrift    = 1
      nkinetic  = 1     
      nfluid    = ncell
      nimpurity = 1
      ALLOCATE(tube    (ntube ))
      ALLOCATE(cell    (ncell ))
      ALLOCATE(field   (nfield))
      ALLOCATE(pin     (npin     ,nion))
      ALLOCATE(photon  (nphoton  ,nion))
      ALLOCATE(drift   (ndrift   ,nion))
      ALLOCATE(kinetic (nkinetic ,nion))
      ALLOCATE(fluid   (nfluid   ,nion))
      ALLOCATE(impurity(nimpurity,nion)) 
c...  Reference plasma solution:
      ref_nion   = 1
      ref_ntube  = 1
      ref_nfluid = 1
      ALLOCATE(ref_tube(ref_ntube))
      ALLOCATE(ref_fluid(ref_nfluid,ref_nion))
c...  Copy DIVIMP grid:
      ion = 1
      ntube = 0
      ncell = 0
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1

        cind1 = ncell + 1
        cind2 = ncell + ike

        ntube = ntube + 1
        tube(ntube)%ir   = ir
        tube(ntube)%type = ringtype(ir)
        tube(ntube)%ikti = ikti2(ir)
        tube(ntube)%ikto = ikto2(ir)
        tube(ntube)%n    = ike
        tube(ntube)%smax = ksmaxs(ir)
        tube(ntube)%pmax = kpmaxs(ir)
        tube(ntube)%rho  = rho(ir,CELL1)
        tube(ntube)%psin = psitarg(ir,1)
        tube(ntube)%cell_index(1) = cind1  ! The "1" must be consistent with "LO" in mod_osm
        tube(ntube)%cell_index(2) = cind2  ! and same for the "2"...
        IF (ir.LT.irsep) THEN
          tube(ntube)%bratio = 0.0
          tube(ntube)%dds    = 0.0
          tube(ntube)%rp     = 0.0
          tube(ntube)%costet = 0.0
        ELSE
          id = idds(ir,2)
          tube(ntube)%bratio(1) = bratio(1,ir)
          tube(ntube)%dds   (1) = dds2(id)
          tube(ntube)%rp    (1) = rp(id)
          tube(ntube)%costet(1) = costet(id)
          tube(ntube)%metric(1) = thetat(id)
          id = idds(ir,1)
          tube(ntube)%bratio(2) = bratio(ike,ir)
          tube(ntube)%dds   (2) = dds2(id)
          tube(ntube)%rp    (2) = rp(id)
          tube(ntube)%costet(2) = costet(id)
          tube(ntube)%metric(2) = thetat(id)
        ENDIF
        
        ncell = cind2
        fluid(cind1:cind2,ion)%te = ktebs(1:ike,ir)
        fluid(cind1:cind2,ion)%ni = knbs (1:ike,ir)
        fluid(cind1:cind2,ion)%vi = kvhs (1:ike,ir)
        fluid(cind1:cind2,ion)%ti = ktibs(1:ike,ir)

        cell(cind1:cind2)%cencar(1) = rs(1:ike,ir)
        cell(cind1:cind2)%cencar(2) = zs(1:ike,ir)
        cell(cind1:cind2)%cencar(3) = 0.0
        cell(cind1:cind2)%vol       = kvols  (1:ike,ir)
        cell(cind1:cind2)%s         = kss   (1:ike,ir)
        cell(cind1:cind2)%p         = kps   (1:ike,ir)
        cell(cind1:cind2)%sbnd(1)   = ksb   (0:ike-1,ir)
        cell(cind1:cind2)%sbnd(2)   = ksb   (1:ike  ,ir)
        cell(cind1:cind2)%pbnd(1)   = kpb   (0:ike-1,ir)
        cell(cind1:cind2)%pbnd(2)   = kpb   (1:ike  ,ir)
        cell(cind1:cind2)%metric    = thetag(1:ike,ir)

        field(cind1:cind2)%bratio = bratio(1:ike,ir)
        IF (load_bfield_data) THEN
          field(cind1:cind2)%b    = bfield_data(1:ike,ir,2)
          field(cind1:cind2)%br   = bfield_data(1:ike,ir,3)
          field(cind1:cind2)%bphi = bfield_data(1:ike,ir,4)
          field(cind1:cind2)%bz   = bfield_data(1:ike,ir,5)
        ELSE
          field(cind1:cind2)%b    = 0.0
          field(cind1:cind2)%br   = 0.0
          field(cind1:cind2)%bphi = 0.0
          field(cind1:cind2)%bz   = 0.0
        ENDIF

        DO ik = 1, ike
          cell(cind1+ik-1)%ik = ik
          cell(cind1+ik-1)%ds = ksb(ik,ir) - ksb(ik-1,ir)
        ENDDO
c        pin(cind1:cind2,ion)%ion = pinion(1:ike,ir)
c        pin(cind1:cind2,ion)%rec = pinrec(1:ike,ir)
c        pin(cind1:cind2,ion)%mom = pinmp (1:ike,ir)
      ENDDO
c...  Set grid quantities:
      grid%n    = ntube
      grid%isep = irsep  - 1
      grid%ipfz = irtrap - 2
      grid%ikto = ikto
      grid%ikti = ikti

c...  Geometry:
      ngrp = 1
      grp(1)%origin = GRP_MAGNETIC_GRID  ! *** NEED AddGroup and CopyGroup functions ***
      grp(1)%type   = GRP_QUADRANGLE
c...  
      nobj = 0
      nsrf = 0
      nvtx = 0
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1

        DO ik = 1, ike
          id = korpg(ik,ir)         
          newobj%group         = ngrp
          newobj%index(IND_IK) = ik
          newobj%index(IND_IR) = ir
          newobj%index(IND_IS) = 0
          newobj%index(IND_CELL   ) = nobj + 1
          newobj%index(IND_FLUID  ) = nobj + 1
          newobj%index(IND_KINETIC) = nobj + 1
          newobj%index(IND_NEUTRAL) = nobj + 1
          newobj%index(IND_FIELD  ) = nobj + 1
          newobj%segment(1) = 0
          newobj%phi        = 0.0
          newobj%nside      = 4
          DO iside = 1, 4
            newsrf%type = SPR_LINE_SEGMENT
            newsrf%obj  = nobj + 1
            newsrf%side = iside
            newsrf%nvtx = 2
            i1 = iside
            i2 = iside + 1
            IF (i2.EQ.5) i2 = 1
            a(1) = DBLE(rvertp(i1,id))
            a(2) = DBLE(zvertp(i1,id))
            a(3) = 0.0D0
            newsrf%ivtx(1) = AddVertex(a) 
            a(1) = DBLE(rvertp(i2,id))
            a(2) = DBLE(zvertp(i2,id))
            a(3) = 0.0D0
            newsrf%ivtx(2) = AddVertex(a) 
            newobj%iside(iside) = AddSurface(newsrf)
          ENDDO
          idum1 = AddObject(newobj)
        ENDDO
      ENDDO

c...  Build connection map:
      CALL BuildConnectionMap(1,nobj)

c...  Setup the old format SOL28 input data:
      osmns28 = osmnnode
      osms28 = 0.0
      DO i1 = 1, osmnnode
        osms28(i1,1) = osmnode(i1)%type
        osms28(i1,5) = osmnode(i1)%rad_x
        osms28(i1,6) = osmnode(i1)%rad_y
      ENDDO

c...  Save legacy data:
      CALL SaveLegacyData('osm_legacy.raw')

c.... Clear arrays:
      IF (ALLOCATED(bfield_data)) DEALLOCATE(bfield_data)

      RETURN
 95   CALL ER('MapRingsToTubes','Error accessing B-field data file',*99)
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE MapTubestoRings(irstart,irend)
      USE mod_sol28_global
      IMPLICIT none
 
      INTEGER irstart,irend,cind1,cind2,ike

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER ir,itube,ion,in

      ion = 1

      DO ir = irstart, irend
        IF (idring(ir).EQ.BOUNDARY) CYCLE

        DO itube = 1, ntube
          IF (tube(itube)%ir.NE.ir) CYCLE

c...      Need to do volume weighted averaging...

          cind1 = tube(itube)%cell_index(1) 
          cind2 = tube(itube)%cell_index(2)

          IF (ir.GE.irsep) THEN
            in = idds(ir,2)
            knds (in) = tube(itube)%ni(1,ion)
            kvds (in) = tube(itube)%vi(1,ion) 
            kteds(in) = tube(itube)%te(1)
            ktids(in) = tube(itube)%ti(1,ion)

            in = idds(ir,1)
            knds (in) = tube(itube)%ni(2,ion)
            kvds (in) = tube(itube)%vi(2,ion)       
            kteds(in) = tube(itube)%te(2)
            ktids(in) = tube(itube)%ti(2,ion)

            ike = nks(ir)
          ELSE
            ike = nks(ir) - 1
          ENDIF

          knes (1:ike,ir) = fluid(cind1:cind2,ion)%ne
          knbs (1:ike,ir) = fluid(cind1:cind2,ion)%ni
          kvhs (1:ike,ir) = fluid(cind1:cind2,ion)%vi
          ktebs(1:ike,ir) = fluid(cind1:cind2,ion)%te
          ktibs(1:ike,ir) = fluid(cind1:cind2,ion)%ti

          osmion(1:ike,ir) = fluid(cind1:cind2,ion)%parion
          osmrec(1:ike,ir) = fluid(cind1:cind2,ion)%parrec
          osmcfp(1:ike,ir) = fluid(cind1:cind2,ion)%parano

          osmmp (1:ike,ir) = fluid(cind1:cind2,ion)%momsrc

          osmcfe(1:ike,ir) = fluid(cind1:cind2,ion)%eneano
          osmqe (1:ike,ir) = fluid(cind1:cind2,ion)%eneion

c...      Finish off core rings:
          IF (ir.LT.irsep) THEN
            knes (nks(ir),ir) = knes (1,ir)
            knbs (nks(ir),ir) = knbs (1,ir)
            kvhs (nks(ir),ir) = kvhs (1,ir) 
            ktebs(nks(ir),ir) = ktebs(1,ir)
            ktibs(nks(ir),ir) = ktibs(1,ir)
          ENDIF

        ENDDO
      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE SetupSOL28Options
      USE mod_sol28_global
      IMPLICIT none


      INCLUDE 'params'
      INCLUDE 'slcom'




      RETURN
 99   STOP
      END

c
c ======================================================================
c
      SUBROUTINE AssignSOL28Nodes_Old(itube,nnode,mnode,node)
      USE mod_sol28
      USE mod_sol28_global
      IMPLICIT none

      INTEGER itube,nnode,mnode       
      TYPE(type_node) :: node(*)

c      SUBROUTINE FindS28Parameters_V4(ir,te,ne,nf,s,new)

      INTEGER ir
      LOGICAL new
      REAL    te(0:6),ne(0:6),nf,s(0:6),isat(0:6)

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      REAL    GetRelaxationFraction,GetJsat

      INTEGER i0,i1,i2,i3,ik,id,index,mode,ik1,ir1,id1,
     .        ikcell(3),ircell(3),id7
      LOGICAL tc,nc,density,tetarget,firsttime
      REAL    frac,t0,t1,n0,n1,A,B,C,tetmp1,tetmp2,coord,expon,
     .        psin0,psin1,prb1,tmp1,val,val0,val1,val2
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd
      DATA firsttime /.FALSE./

c...  
      ir = tube(itube)%ir
c      IF (output) WRITE(0,*) ' ASSIGN PARAMS: IR=',ir
     

      s (0) = 0.0
      ne(1) = 0.0
      te(1) = 0.0
      te(2) = 0.0
      te(3) = 0.0
      te(4) = 0.0
      te(5) = 0.0


      frac = GetRelaxationFraction()


      IF (new.AND.frac.NE.0.0) THEN
c        WRITE(0,*) 'FRAC:',frac,rel_step,rel_nstep
      ENDIF


      DO i1 = 2, osmns28
        i0 = i1 - 1
        IF ((osms28(i0,1).NE.osms28(i1,1)).OR.osms28(i1,2).EQ.0.0) CYCLE

c...    Do not apply data if IR is outside specified range:
        IF ((osms28(i1,11).NE.0.0.AND.REAL(ir).LT.osms28(i1,11) ).OR.
     .      (osms28(i1,12).NE.0.0.AND.REAL(ir).GT.osms28(i1,12)))
     .    CYCLE

c...    Check that rings from different grid regions are not in the same group
c       of rings:
        DO i2 = NINT(osms28(i1,11)), NINT(osms28(i1,12))-1
          IF (ringtype(i2).NE.ringtype(i2+1)) THEN
c            WRITE(0,*) 'RINGTYPE:',i1
c            WRITE(0,*) 'RINGTYPE:',i2,ringtype(i2)
c            WRITE(0,*) 'RINGTYPE:',i2+1,ringtype(i2+1)
c            WRITE(0,*) 'RINGTYPE:',osmns28
c            CALL ER('FindS28Parameters_V3','Mixed RINGTYPE',*99)
            IF (sloutput.AND.firsttime) THEN
              firsttime = .FALSE.
              WRITE(0,*)
              WRITE(0,*) '-----------------------------------------'
              WRITE(0,*) ' THAT FUNNY THING ABOUT MIXING REGIONS!? '
              WRITE(0,*) '-----------------------------------------'
              WRITE(0,*)
            ENDIF
          ENDIF      
        ENDDO


        IF (.FALSE..AND.osms28(i1,2).EQ.99.0) THEN
        ELSE

c MODE                          P1           P2
c   1 - power law               coordinate   index
c   2 - exponential v0-v2       coordinate   index
c   3 - exponential to infinity
c   4 - from probe data         coordinate   probe number
c
c
c   coord = 1 - linear on line segment
c         = 2 - linear on line segment, but from first to last ring intersection
c         = 3 - PSIn 
c         = 4 - RHO
c

          index = NINT(osms28(i1,1)) - 1
          mode  = NINT(osms28(i1,2))
          coord = osms28(i1,3)
          expon = osms28(i1,4)

c...  Decide if specified upstream data is density or pressure:
          density = .TRUE.   ! *** I DON'T LIKE THIS SOLUTION ***
          IF (index.EQ.3.AND.
     .        (ringtype(ir).EQ.SOL1.AND.s28nemode   .EQ.1.OR.
     .         ringtype(ir).EQ.PFZ .AND.s28nemodepfz.EQ.1)) 
     .      density = .FALSE.

          a1 = DBLE(osms28(i1-1,5))
          a2 = DBLE(osms28(i1-1,6))
          b1 = DBLE(osms28(i1  ,5))
          b2 = DBLE(osms28(i1  ,6))

          DO ik = 1, nks(ir)

            id = korpg(ik,ir)
            c1 = 0.5D0 * DBLE(rvertp(1,id) + rvertp(2,id))
            c2 = 0.5D0 * DBLE(zvertp(1,id) + zvertp(2,id))
            d1 = 0.5D0 * DBLE(rvertp(3,id) + rvertp(4,id))
            d2 = 0.5D0 * DBLE(zvertp(3,id) + zvertp(4,id))
       
            CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
 
            IF (tab.GE.0.0.AND.tab.LT.1.0.AND.
     .          tcd.GE.0.0.AND.tcd.LT.1.0) THEN
c...          Intersecting between the line segment and the ring is found:

c...          These are here in case psin0 gets picked up when linking exponential
c             decay data to neighbouring ring:
              IF (ringtype(ir).EQ.PFZ) THEN
                psin0 = -HI
                psin1 =  HI
              ELSE
                psin0 = HI
                psin1 = HI
              ENDIF

              IF (.FALSE.) THEN
c...            Need a check to know if the ring is over-constrained:
              ELSEIF (index.GE.0.AND.index.LE.6) THEN
c              ELSEIF (index.GE.1.AND.index.LE.5) THEN

                IF (index.GE.1.AND.index.LE.5) THEN
                  s(index)=ksb(ik-1,ir)+
     .                     SNGL(tcd)*(ksb(ik,ir)-ksb(ik-1,ir))
                ENDIF

c...            Find data boundary values -- NEEDS WORK!:
                i2 = i0

                i3 = i1

c...            Flag...
                tetarget = .FALSE.
                IF ((index.EQ.1.OR.index.EQ.5).AND.osms28(i3,7).LT.0.0) 
     .            tetarget = .TRUE.

                IF     (mode.EQ.1.OR.mode.EQ.2.OR.mode.EQ.3) THEN
c...              Interpolation boundary values provided in the input file:
c                    WRITE(0,*) 'DATA1:',index,osms28(i2,7),osms28(i2,8)

                  t0 = osms28(i2,7) + frac*(osms28(i2,9) -osms28(i2,7))
                  n0 = osms28(i2,8) + frac*(osms28(i2,10)-osms28(i2,8))
                  t1 = osms28(i3,7) + frac*(osms28(i3,9) -osms28(i3,7))
                  n1 = osms28(i3,8) + frac*(osms28(i3,10)-osms28(i3,8))

                  IF ((index.EQ.3.OR.index.EQ.4.OR.index.EQ.5).AND.
     .                (osms28(i2,7).EQ.-99.0.OR.
     .                 osms28(i2,8).EQ.-99.0)) THEN
c...                Linking to another plasma region that is *already* calculated: 
                    CALL FindCell(i2,i3,ir,ikcell,ircell)
                    IF (ringtype(ir).EQ.PFZ) THEN
                      ik1 = ikouts(ikcell(2),ircell(2))
                      ir1 = irouts(ikcell(2),ircell(2))
c                      WRITE(0,*) 'FOUND YOU PFZ:',ik1,ir1,
c     .  density,tetarget,index
                    ELSE
                      ik1 = ikins(ikcell(2),ircell(2))
                      ir1 = irins(ikcell(2),ircell(2))
                    ENDIF
                    IF (osms28(i2,7).EQ.-99.0) t0 = ktebs(ik1,ir1)
                    IF (osms28(i2,8).EQ.-99.0) THEN
                      IF (density) THEN
                        n0 = knbs (ik1,ir1)
c                        WRITE(0,*) 'denisity POWER:',n0
                      ELSE
                        n0 = 2.0 * ktebs(ik1,ir1) * knbs(ik1,ir1)  ! Add Ti and M? 
                      ENDIF
                    ENDIF
c...                Shouldn't really be outer target (all this would go away if PSITARG was
c                   assigned properly):
                    IF (coord.EQ.3) psin0 = psitarg(ir1,1)
c                   IF (coord.EQ.4) rho0 =
c                    WRITE(0,*) 'DATA2:',t0,psin0,ik1,ir1,index
c                    WRITE(0,*) 'DATA2:',n0,psin0,ik1,ir1
c                    STOP 'rgsd'
                  ENDIF

c...              Make sure that t0,1 are positive:
                  IF (tetarget) THEN
                    t0 = ABS(t0)
                    t1 = ABS(t1)
                  ENDIF

                  IF (coord.EQ.1) THEN
c...                Linear along the line segment, nice and simple:
                    val0 = 0.0
                    val1 = 1.0
                    val = SNGL(tab)
c                    val = SNGL(tcd) ! BUG
                  ELSE
c...              
                    IF (NINT(osms28(i2,11)).NE.NINT(osms28(i3,11)).OR.
     .                  NINT(osms28(i2,12)).NE.NINT(osms28(i3,12)).OR.
     .                  NINT(osms28(i2,11)).GT.NINT(osms28(i2,12)).OR.
c...                    This is here because PSITARG is currently not assigned in the core... 
     .                  NINT(osms28(i2,11)).LT.irsep)
     .                CALL ER('FindS28Parameters_V3','Invalid '//
     .                        'IR range',*99)

c...                Need range of PSIn over the segment:
                    DO ir1 = NINT(osms28(i2,11)), NINT(osms28(i2,12)) 
                      DO ik1 = 1, nks(ir1)
                        id1 = korpg(ik1,ir1)
                        c1 = 0.5D0 * DBLE(rvertp(1,id1) + rvertp(2,id1))
                        c2 = 0.5D0 * DBLE(zvertp(1,id1) + zvertp(2,id1))
                        d1 = 0.5D0 * DBLE(rvertp(3,id1) + rvertp(4,id1))
                        d2 = 0.5D0 * DBLE(zvertp(3,id1) + zvertp(4,id1))
                        CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
                        IF (tab.GE.0.0.AND.tab.LT.1.0.AND.
     .                      tcd.GE.0.0.AND.tcd.LT.1.0) THEN
c...                      Should be more careful which PSITARG is used and when...
c                          WRITE(0,*) 'MARK:',ir1,psin0,psin1
                          IF     (ringtype(ir).EQ.SOL1) THEN
                            IF (psin0.EQ.HI) psin0 = psitarg(ir1,1)
                            IF (psin0.NE.HI) psin1 = psitarg(ir1,1)
c                            WRITE(0,*) 'PSIN:',ir1,psin0,psin1
                          ELSEIF (ringtype(ir).EQ.PFZ) THEN
                            psin0 = MAX(psin0,psitarg(ir1,1))
                            psin1 = MIN(psin1,psitarg(ir1,1))
                          ELSE
                            CALL ER('FindS28Parameters_V3','Invalid '//
     .                              'RINGTYPE',*99)
                          ENDIF
                        ENDIF
                      ENDDO
                    ENDDO


c                    WRITE(0,*) 'PSIn:',ir,psin0,psin1

                  
                    IF     (coord.EQ.2) THEN
c...                  Spatial along the IR-range of applicability:
                      STOP 'NOT READY 2'
                    ELSEIF (coord.EQ.3) THEN
c...                  PSIn:
                      val0 = 0.0
                      val1 = ABS(psin1         - psin0)
                      val  = ABS(psitarg(ir,1) - psin0)
                    ELSEIF (coord.EQ.4) THEN
c...                  RHO:
                      STOP 'NOT READY 4'
                    ELSE
                      WRITE(0,*) 'I0,I1=',i0,i1
                      CALL ER('S28params_v3','Invalid COORD A',*99)
                    ENDIF
                  ENDIF

                ELSEIF (mode.EQ.4) THEN
c...              Load probe data, dummy values here:
                  t0 = osms28(i2,7)
                  t1 = t0
                  n0 = osms28(i2,8)
                  n1 = n0

                ELSE
                  CALL ER('S28params_v3','Invalid MODE',*99)   
                ENDIF

c...            Check if quantities should be assigned:
                tc = .TRUE.
                nc = .TRUE.
                IF (t0.EQ.0.0.OR.t1.EQ.0.0) tc = .FALSE.
                IF (n0.EQ.0.0.OR.n1.EQ.0.0) nc = .FALSE.

c                IF (ir.GT.22.AND.ir.LT.irwall.AND.index.EQ.1.AND.
c     .              .NOT..FALSE..AND.(.TRUE..OR.index.EQ.4)) THEN
c                  WRITE(0,*) 'PSIn:',ir,psitarg(ir,1),psin0,psin1
c                  WRITE(0,*) 'i2,3:',i2,i3
c                  WRITE(0,*) 'val :',ir,val,val0,val1
c                  WRITE(0,*) 'Te  :',ir,t0,t1
c                  WRITE(0,*) 'ne  :',ir,n0,n1
c                ENDIF

                IF     (mode.EQ.1) THEN
c...              Power law between v1 and v2:
                  val2 = (val - val0) / (val1 - val0)
c                  IF (index.EQ.5) WRITE(0,*) '5:',val2,ir

                  IF (tc) te(index) = t0 + val2**expon * (t1 - t0)
                  IF (nc) ne(index) = n0 + val2**expon * (n1 - n0)

c                  IF (index.EQ.3.AND.ir.EQ.71) THEN
c                    WRITE(0,*) 'A:',val,n0,n1,density
c                  ENDIF
      

                ELSEIF (mode.EQ.2) THEN
c...              Exponential decay between v1 and v2:
                  C = expon  ! BUG!
c                  C = -expon
c                  val = 1.0

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



                ELSEIF (mode.EQ.3) THEN
c...              Exponential decay to infinity:
                  C = expon
                  A = t0 - t1
                  B = t1 
                  IF (tc) te(index) = A * EXP(-val / C) + B
                  A = n0 - n1
                  B = n1
                  IF (nc) ne(index) = A * EXP(-val / C) + B

                ELSEIF (mode.EQ.4) THEN
c...              Load probe data from the .experiments file:
                  IF     (coord.EQ.3) THEN
                    prb1 = -1.0
                  ELSEIF (coord.EQ.4) THEN
                    prb1 = -3.0
                  ELSE
                      WRITE(0,*) 'I0,I1=',i0,i1,coord
                    CALL ER('S28params_v3','Invalid COORD B',*99)   
                  ENDIF

                  IF (tc) THEN
                    tmp1 = prb1
                    CALL LoadProbeDataS28(ir,NINT(osms28(i2,7)),2,tmp1) 
                    te(index) = tmp1
                    WRITE(0,*) 'TE PROBE:',tmp1,ir,psitarg(ir,1)
                    IF     (expon.EQ.1.0) THEN
                      STOP 'OPTION NOT READY'
                    ELSEIF (expon.EQ.2.0) THEN
c...                  
                      te(0) = -te(index)
                      te(6) = -te(index)
                    ELSEIF (expon.EQ.3.0) THEN
c...                  
                      te(0) = 97.0
                      te(6) = 97.0
                    ENDIF
                  ENDIF

                  IF (nc) THEN
                    tmp1 = prb1
                    CALL LoadProbeDataS28(ir,NINT(osms28(i2,8)),1,tmp1) 
                    ne(index) = tmp1
                    WRITE(0,*) 'NE PROBE:',tmp1,ir,psitarg(ir,1)
                    IF     (expon.EQ.1.0) THEN
                      STOP 'OPTION NOT READY'
                    ELSEIF (expon.EQ.2.0) THEN
c...                  
                      ne(0) = -ne(index)
                      ne(6) = -ne(index)
                    ELSEIF (expon.EQ.3.0) THEN
                    ENDIF
                  ENDIF

c                  WRITE(0,*) 'PROBIN:',prb1,te(index),ne(index),tc,nc

                ELSE
                  CALL ER('S28params_v3','Invalid MODE',*99)   
                ENDIF

              ELSE
                CALL ER('FindSOL28Parameters','Invalid parameter '//
     .                  'index',*99)
              ENDIF


              IF (tetarget) THEN
c                WRITE(0,*) 'TETARGET:',tetarget,te(index),i3,index
                te(index) = -te(index)
              ENDIF

              IF (.NOT.density.AND.index.EQ.3.AND.nc.AND.mode.NE.4) THEN
c...            Convert the pressure value to density:
                IF (te(3).EQ.0.0) THEN
                  CALL ER('FindSOL28Parameters','Te3 not assigned',*99)
                ELSE
c...              Assumes that the Mach no. is low (note: Ti.NE.Te is not a problem
c                 since that situation is currently resolved in the SOL28 routine):
c                  WRITE(0,*) 'CONVERTING PRESSURE TO DENSITY',ne(3)
c                  WRITE(0,*) 'OLD D:',ne(3)
                  ne(3) = ne(3) / (2.0 * te(3))
c                  WRITE(0,*) 'NEW D:',ne(3)
                ENDIF
              ENDIF

            ENDIF

          ENDDO

c          IF (.NOT.density.AND.index.EQ.3.AND.nc.AND.mode.NE.4) THEN
cc...        Convert the pressure value to density:
c            IF (te(3).EQ.0.0) THEN
c              CALL ER('FindSOL28Parameters','Te(3) not assigned',*99)
c            ELSE
cc...          Assumes that the Mach no. is low (note: Ti.NE.Te is not a problem
cc             since that situation is currently resolved in the SOL28 routine):
cc              WRITE(0,*) 'CONVERTING PRESSURE TO DENSITY',ne(3)
c              WRITE(0,*) 'OLD D:',ne(3)
c              ne(3) = ne(3) / (2.0 * te(3))
c              WRITE(0,*) 'NEW D:',ne(3)
c            ENDIF
c          ENDIF

c...    End of OSMS28(I1,1).EQ.98.0 block:      
        ENDIF

      ENDDO




c          WRITE(0,*) 'NE3=',ne(3)
      IF (.TRUE.) THEN
c...    Specify te(1) from target data.  The -ve is to trigger
c       te(0)=te(1):
        IF (te(1).EQ.-98.0) te(1) = -kteds(idds(ir,2))

c...    Specify te(5) from target data.  The -ve is to trigger
c       te(6)=te(5):
        IF (te(5).EQ.-98.0) te(5) = -kteds(idds(ir,1))
      ENDIF 


     

      IF (.NOT..TRUE.) THEN


        IF (te(3).EQ.-97.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,2,tetmp1) 
          te(3) = tetmp1
c...      Don't necessarily want both targets done...:
          te(0) = 97.0
          te(6) = 97.0
        ENDIF
        IF (ne(3).EQ.-97.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,1,tetmp1) 
          ne(3) = tetmp1
        ENDIF


        IF (te(3).EQ.-98.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,2,tetmp1) 
c          CALL LoadProbeDataS28(ir,77,2,tetmp2) 
c          te(3) = 0.50 * tetmp1 + 0.50 * tetmp2
          te(3) = tetmp1
          te(0) = -te(3)
          te(6) = -te(3)
        ENDIF
        IF (ne(3).EQ.-98.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,1,tetmp1) 
c          CALL LoadProbeDataS28(ir,77,1,tetmp2) 
c          ne(3) = 0.50 * tetmp1 + 0.50 * tetmp2
          ne(3) = tetmp1
          ne(0) = -ne(3)
          ne(6) = -ne(3)
        ENDIF

c        WRITE(0,*) 'NE(3):',ne(3),te(3)

        IF (te(3).EQ.-99.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,70,2,tetmp1) 
          CALL LoadProbeDataS28(ir,77,2,tetmp2) 
c         te(3) = 1.00 * tetmp1 + 0.00 * tetmp2
          te(3) = 0.50 * tetmp1 + 0.50 * tetmp2
        ENDIF
        IF (ne(3).EQ.-99.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,70,1,tetmp1) 
          CALL LoadProbeDataS28(ir,77,1,tetmp2) 
c          ne(3) = 1.00 * tetmp1 + 0.00 * tetmp2
          ne(3) = 0.50 * tetmp1 + 0.50 * tetmp2
          WRITE(PINOUT,*) 'COMPROMISE:',ir,psitarg(ir,2),ne(3),te(3)
        ENDIF

        IF (te(3).LE.-1.0) CALL LoadProbeDataS28(ir,s28probe,2,te(3))
c *HACKISH*
        IF (ne(3).LE.-1.0.AND.ne(3).GT.-100.0) 
     .    CALL LoadProbeDataS28(ir,s28probe,1,ne(3))

        IF (ne(3).LT.-100.0) THEN
c...      Pressure specified, not density.  Calculate ne(3) from the pressure
c         and the specified te(3) value, assuming Ti=Te:
          ne(3) = -ne(3) / (2.0 * te(3))
        ENDIF

c...    Specify te(1) from target data.  The -ve is to trigger
c       te(0)=te(1):
        IF (te(1).EQ.-99.0) te(1) = -kteds(idds(ir,2))

c...    Specify te(5) from target data.  The -ve is to trigger
c       te(6)=te(5):
        IF (te(5).EQ.-99.0) te(5) = -kteds(idds(ir,1))

      ELSE
      ENDIF


      IF (.TRUE.) THEN

        s(0) = 0.0
        s(6) = ksmaxs(ir)

        DO i1 = 0, 6
          DO ik = 1, nks(ir)
            IF (s(i1).GE.ksb(ik-1,ir).AND.s(i1).LE.ksb(ik,ir)) THEN
              node(i1+1)%icell = ik
              IF (i1.EQ.3) s(3) = kss(ik,ir)
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

        id1 = idds(ir,2)
        id7 = idds(ir,1)
        IF (node(1)%ne.EQ.0.0) node(1)%ne = knds (id1)
        IF (node(7)%ne.EQ.0.0) node(7)%ne = knds (id7)
        IF (node(1)%te.EQ.0.0) node(1)%te = kteds(id1)
        IF (node(7)%te.EQ.0.0) node(7)%te = kteds(id7)
        IF (node(1)%ti(1).EQ.0.0) node(1)%ti(1) = ktids(id1)
        IF (node(7)%ti(1).EQ.0.0) node(7)%ti(1) = ktids(id7)

c REAL FUNCTION GetJsat(te,ti,ne,v)
        id = id1
        node(1)%jsat(1)= GetJsat(kteds(id),ktids(id),knds(id),kvds(id))
        node(1)%jsat(1) = -ABS(node(1)%jsat(1))
        id = id7
        node(7)%jsat(1)= GetJsat(kteds(id),ktids(id),knds(id),kvds(id))
        node(7)%jsat(1) = -ABS(node(7)%jsat(1))


      ENDIF


      RETURN
99    CONTINUE
      WRITE(0,*) 'IK,IR=',ik,ir,i1,osms28(i1,1)
      WRITE(0,*) 'TAB,TCD=',tab,tcd
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE SaveLegacyData(fname)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      CHARACTER*(*) fname
      INTEGER fp,i1,i2

      fp = 99
      OPEN(UNIT=fp,FILE=fname(1:LEN_TRIM(fname)),ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='REPLACE',ERR=98)            
      WRITE(fp,ERR=98) 1.00 
      WRITE(fp,ERR=98) MAXNRS
      WRITE(fp,ERR=98) (grdntreg(i1),grdntseg(1:MAXNRS,i1),
     .                  (grdtseg(1:MAXNRS,i2,i1),i2=1,MAXNRS),i1=1,2)

      CLOSE (fp)
      
      RETURN
 98   CALL ER('SaveGrid','Problem saving legacy data file',*99)
 99   STOP
      END
c
c ======================================================================
c
c      SUBROUTINE SetTargetConditions_Legacy
c...  Dummy
c      END
c
c ======================================================================
c


      SUBROUTINE BuildLinearGrid
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'


      INTEGER id,ik,ir,i1,grid_option
      REAL*4  r,delr,L,r1,r2,z1,z2,delta,frac,frac1,frac2,vessel_radius


      grid_option = 6

      vessel_radius = 0.02  
c      vessel_radius = 0.12  ! ULS


      SELECTCASE (grid_option)
        CASE (1)  ! Full vessel, mirrored
          L = 3.6                     ! Total length of mirrored plasma column (m)
          r = 0.015                   ! Plasma radius (m)
          delr = (vessel_radius - r)  ! Distance from plasma to outer wall (m)
          maxrings = 10               ! Number of flux tubes (if changed, also need to change triangle grid in input file)
          nks(1:maxrings) = 100       ! Number of cells on each tube
        CASE (2)  ! Full vessel
          L = 1.8
          r = 0.015          
          delr = (vessel_radius - r)  
          maxrings = 10      
          nks(1:maxrings) = 100         
        CASE (3) ! Target chamber
          L = 0.56
          r = 0.015          
          delr = (vessel_radius - r)  
          maxrings = 10      
          nks(1:maxrings) = 50  ! 175
        CASE (4) ! Target chamber: fancy
          vessel_radius = 0.16 
          L = 0.55
          r = 0.08          
          delr = (vessel_radius - r)  
          maxrings = 50      
          nks(1:maxrings) = 150
        CASE (5) ! Target chamber: fancy #2, full vessel and small volume 
          vessel_radius = 0.16 
          L = 0.55
          r = 0.03          
          delr = (vessel_radius - r)  
          maxrings = 50      
          nks(1:maxrings) = 150
        CASE (6) ! Target chamber: fancy #3, small volume 
          vessel_radius = 0.05 
          L = 0.55
          r = 0.03          
          delr = (vessel_radius - r)  
          maxrings = 50      
          nks(1:maxrings) = 150
      ENDSELECT

      id = 0

      r0 = 0.0000001  ! Need this tiny displacement to keep EIRENE04 from falling over 
      z0 = 0.0

      DO ir = 1, maxrings

        IF (.TRUE.) THEN
          IF (grid_option.EQ.4) THEN
            frac1 = (REAL(ir-1) / REAL(maxrings))**0.7
            frac2 = (REAL(ir  ) / REAL(maxrings))**0.7
            r1 = r0 + frac1 * r
            r2 = r0 + frac2 * r
          ELSE
            frac  = REAL(ir-1) / REAL(maxrings)
            delta = r / REAL(maxrings)
            r1 = r0 + frac * r 
            r2 = r1 + delta
          ENDIF
        ENDIF

C       Krieger IPP/07 - SUN compiler does not know SNGL, replaced by REAL  -strange since SNGL is used elsewhere... -SL
C       psitarg(ir,1) = ABS(0.5*(SNGL(r1+r2)))
C       psitarg(ir,2) = ABS(0.5*(SNGL(r1+r2)))
        psitarg(ir,1) = ABS(0.5*(REAL(r1+r2)))
        psitarg(ir,2) = ABS(0.5*(REAL(r1+r2)))
        idring(ir) = TARTOTAR

c        WRITE(0,*) 'IR:',ir,psitarg(ir,1)

c       nks(ir) = 100

        DO ik = 1, nks(ir)

          SELECTCASE (grid_option)
            CASE (1)  ! Full vessel, mirrored
              IF (.TRUE.) THEN
                frac = REAL(ik-1) / REAL(nks(ir)) 
                delta = L / REAL(nks(ir)) 
                z1 = (0.5 - frac) * L
                z2 = z1 - delta
              ENDIF
            CASE (2)  ! Full vessel
              IF (.TRUE.) THEN
                frac = REAL(ik-1) / REAL(nks(ir)) 
                delta = L / REAL(nks(ir)) 
                z1 = (1.0 - frac) * L
                z2 = z1 - delta       
              ENDIF
            CASE (3)  ! Target chamber
              IF (.TRUE.) THEN
                frac = REAL(ik-1) / REAL(nks(ir)) 
                delta = L / REAL(nks(ir)) 
c                z1 = frac * L
c                z2 = z1 + delta       
                z1 = (1.0 - frac) * L
                z2 = z1 - delta       
              ENDIF
            CASE (4:6) ! Target chamber: fancy
              frac1 = REAL(ik-1) / REAL(nks(ir)) 
              frac2 = REAL(ik  ) / REAL(nks(ir)) 
              frac1 = SIGN(0.5,frac1-0.5)*(ABS(frac1-0.5)/0.5)**1.00+0.5
              frac2 = SIGN(0.5,frac2-0.5)*(ABS(frac2-0.5)/0.5)**1.00+0.5
              z1 = (1.0 - frac1) * L
              z2 = (1.0 - frac2) * L     
          ENDSELECT

          bratio(ik,ir) = 1.0
          kbfs  (ik,ir) = 1.0
          bts   (ik,ir) = cbphi 

          id = id + 1

          korpg(ik,ir) = id

          nvertp(id) = 4

          IF (ik.EQ.1) THEN
            rvertp(1,id) = r1
            rvertp(2,id) = r2
            zvertp(1,id) = z1
            zvertp(2,id) = z1
          ELSE
            rvertp(1,id) = rvertp(4,id-1)
            rvertp(2,id) = rvertp(3,id-1)
            zvertp(1,id) = zvertp(4,id-1)
            zvertp(2,id) = zvertp(3,id-1)
          ENDIF

          rvertp(3,id) = r2
          rvertp(4,id) = r1
          zvertp(3,id) = z2
          zvertp(4,id) = z2

          rs(ik,ir) = 0.0
          zs(ik,ir) = 0.0
          DO i1 = 1, nvertp(id)
            rs(ik,ir) = rs(ik,ir) + rvertp(1,id)
            zs(ik,ir) = zs(ik,ir) + zvertp(1,id)
          ENDDO
          rs(ik,ir) = rs(ik,ir) / REAL(nvertp(id))
          zs(ik,ir) = zs(ik,ir) / REAL(nvertp(id))

        ENDDO

      ENDDO

      npolyp  = id
      vpolmin = (MAXNKS*MAXNRS - npolyp) / 2 + npolyp
      vpolyp  = vpolmin

      ikto = 2
      ikti = 3

      rves = 0
      rvesm = 0

      irsep  = 1
      irwall = maxrings 
      irtrap = irwall
      nrs    = irwall
      nbr    = 0

      WRITE(0,*) 'NVERT:',nvertp(5)

      CALL InsertRing(1         ,BEFORE,PERMANENT)
      CALL InsertRing(maxrings+1,AFTER ,PERMANENT)

      WRITE(0,*) 'NVERT:',nvertp(5)

c...  Necessary..? 
      cutring = 1
      cutpt1 = ikto
      cutpt2 = ikti

      idring(1) = -1
      idring(nrs) = -1


c...  Modify the grid based on entries in the GRDMOD array assigned 
c     from the input file:
c      IF (grdnmod.NE.0) CALL TailorGrid


      rxp =  0.0 
      zxp = -0.25 * L

      rmin = HI
      rmax = LO
      zmin = HI
      zmax = LO
      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          rmin = MIN(rmin,rs(ik,ir))
          rmax = MAX(rmax,rs(ik,ir))
          zmin = MIN(zmin,zs(ik,ir))
          zmax = MAX(zmax,zs(ik,ir))
        ENDDO
      ENDDO

c...  Neutral wall
      IF (cneur.EQ.4) THEN

         SELECTCASE (grid_option)
           CASE (1)  ! Full vessel, mirrored
             nves = 20
             ir = irwall-1
             r1 = rvertp(2,korpg(1      ,ir))
             r2 = r1 + delr
             z1 = zvertp(2,korpg(1      ,ir))
             z2 = zvertp(3,korpg(nks(ir),ir))
  
             rves(1)  =  r1
             zves(1)  =  z1
             rves(2)  =  r2
             zves(2)  =  z1
  
             rves(3)  =  r2
             zves(3)  =  L / 2.0 - 0.56
             rves(4)  =  r1 + 0.0101 
             zves(4)  =  L / 2.0 - 0.56
             rves(5)  =  r1 + 0.0001 
             zves(5)  =  L / 2.0 - 0.57
             rves(6)  =  r2
             zves(6)  =  L / 2.0 - 0.57
  
             rves(7)  =  r2
             zves(7)  =  0.03
             rves(8)  =  r1 + 0.0001
             zves(8)  =  0.03
             rves(9)  =  r1 + 0.0001
             zves(9)  =  0.02
             rves(10) =  r2
             zves(10) =  0.02
  
             rves(11) =  r2
             zves(11) = -0.02
             rves(12) =  r1 + 0.0001
             zves(12) = -0.02
             rves(13) =  r1 + 0.0001
             zves(13) = -0.03
             rves(14) =  r2
             zves(14) = -0.03
  
             rves(15) =  r2
             zves(15) = -L / 2.0 + 0.56
             rves(16) =  r1 + 0.0001
             zves(16) = -L / 2.0 + 0.56
             rves(17) =  r1 + 0.0101
             zves(17) = -L / 2.0 + 0.57
             rves(18) =  r2
             zves(18) = -L / 2.0 + 0.57
  
             rves(19) =  r2
             zves(19) =  z2
             rves(20) =  r1
             zves(20) =  z2
           CASE (2)  ! Full vessel
             nves = 12
             ir = irwall-1
             r1 = rvertp(2,korpg(1      ,ir))
             r2 = r1 + delr
             z1 = zvertp(2,korpg(1      ,ir))
             z2 = zvertp(3,korpg(nks(ir),ir))
  
             rves(1)  =  r1
             zves(1)  =  z1
             rves(2)  =  r2
             zves(2)  =  z1
  
             rves(3)  =  r2
             zves(3)  =  L - 0.56
             rves(4)  =  r1 + 0.0101 
             zves(4)  =  L - 0.56
             rves(5)  =  r1 + 0.0001 
             zves(5)  =  L - 0.57
             rves(6)  =  r2
             zves(6)  =  L - 0.57
  
             rves(7)  =  r2
             zves(7)  =  0.03
             rves(8)  =  r1 + 0.0001
             zves(8)  =  0.03
             rves(9)  =  r1 + 0.0001
             zves(9)  =  0.02
             rves(10) =  r2
             zves(10) =  0.02
  
             rves(11) =  r2
             zves(11) =  z2
             rves(12) =  r1
             zves(12) =  z2
           CASE (3)  ! Target chamber
             nves = 4
             ir = irwall-1
             r1 = rvertp(2,korpg(1      ,ir))
             r2 = r1 + delr
             z1 = zvertp(2,korpg(1      ,ir))
             z2 = zvertp(3,korpg(nks(ir),ir))
  
             rves(1) =  r1
             zves(1) =  z1
             rves(2) =  r2
             zves(2) =  z1
  
             rves(3) =  r2
             zves(3) =  z2 
             rves(4) =  r1
             zves(4) =  z2
           CASE (4:5)  ! Target chamber: fancy
             nves = 10
             ir = irwall-1
             r1 = rvertp(2,korpg(1      ,ir))
             r2 = r1 + delr
             z1 = zvertp(2,korpg(1      ,ir))
             z2 = zvertp(3,korpg(nks(ir),ir))
  
             rves(1) =  r1
             zves(1) =  z1
             rves(2) =  r2
             zves(2) =  z1

             rves(3) =  r2
             zves(3) =  0.11

             rves(4) =  0.21
             zves(4) =  0.11

             rves(5) =  0.21
             zves(5) =  0.06

             rves(6) =  0.212
             zves(6) =  0.06

             rves(7) =  0.212
             zves(7) =  0.05

             rves(8) =  0.21
             zves(8) =  0.05

             rves(9) =  0.21
             zves(9) =  z2
  
             rves(10) =  r1
             zves(10) =  z2
           CASE (6)  ! Target chamber: fancy #3, small volume
             nves = 10
             ir = irwall-1
             r1 = rvertp(2,korpg(1      ,ir))
             r2 = r1 + delr
             z1 = zvertp(2,korpg(1      ,ir))
             z2 = zvertp(3,korpg(nks(ir),ir))
  
             rves(1) =  r1
             zves(1) =  z1
             rves(2) =  r2
             zves(2) =  z1

             rves(3) =  r2
             zves(3) =  0.11

             rves(4) =  r2 + 0.01
             zves(4) =  0.11

             rves(5) =  r2 + 0.01
             zves(5) =  0.06

             rves(6) =  r2 + 0.015
             zves(6) =  0.06

             rves(7) =  r2 + 0.015
             zves(7) =  0.05

             rves(8) =  r2 + 0.01
             zves(8) =  0.05

             rves(9) =  r2 + 0.01
             zves(9) =  z2
  
             rves(10) =  r1
             zves(10) =  z2
        ENDSELECT
      ENDIF


      IF (.TRUE.) THEN
        nvesm = nves - 1
        DO i1 = 1, nves-1
          rvesm(i1,1) = rves(i1)
          zvesm(i1,1) = zves(i1)
          rvesm(i1,2) = rves(i1+1)
          zvesm(i1,2) = zves(i1+1)
        ENDDO
      ENDIF
 
c      CALL DumpGrid('BUILDING LINEAR GRID')


      CALL OutputData(85,'Linear')


c...  Add virtual boundary cells, which will be stripped off later:
      IF (CTARGOPT.EQ.0.OR.CTARGOPT.EQ.1.OR.CTARGOPT.EQ.2.OR.
     .    CTARGOPT.EQ.3.OR.CTARGOPT.EQ.6) 
     .   CALL AddPoloidalBoundaryCells

c      STOP 'WHA-WHO!'

      RETURN
 99   STOP
      END
