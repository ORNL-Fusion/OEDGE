c     -*-Fortran-*-
c
c ======================================================================
c
      REAL FUNCTION ParticleSource(itube)
      USE mod_geometry
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      REAL CalcFlux

      INTEGER, INTENT(IN) :: itube

      INTEGER cind1,cind2,ic,ip,ion
      REAL    sumtot,sumion,sumrec,sumano,sumflx

      ParticleSource = 0.0

      ion = 1

      cind1 = tube(itube)%cell_index(LO)
      cind2 = tube(itube)%cell_index(HI)
      sumion = 0.0
      DO ic = cind1, cind2
        ip = obj(cell2obj(ic))%index(IND_FLUID)
        sumion = sumion + cell(ic)%vol * 
     .           fluid(ip,ion)%parion
c     .           (fluid(ip,ion)%parion - fluid(ip,ion)%parrec +
c     .            fluid(ip,ion)%parano)
c        IF (itube.EQ.4) 
c     .    WRITE(0,*) 'CHECK:',ic,ip,cell(ic)%vol,
c     .               fluid(ip,ion)%parion,fluid(ip,ion)%parrec
      ENDDO

      IF (itube.EQ.4) 
     .  WRITE(0,*) 'FLUX:',itube,sumion

      sumrec = 0.0
      DO ic = cind1, cind2
        ip = obj(cell2obj(ic))%index(IND_FLUID)
        sumrec = sumrec + cell(ic)%vol * 
     .           fluid(ip,ion)%parrec
c     .           (fluid(ip,ion)%parion - fluid(ip,ion)%parrec +
c     .            fluid(ip,ion)%parano)
c        IF (itube.EQ.4) 
c     .    WRITE(0,*) 'CHECK:',ic,ip,cell(ic)%vol,
c     .               fluid(ip,ion)%parion,fluid(ip,ion)%parrec
      ENDDO

      IF (itube.EQ.4) 
     .  WRITE(0,*) 'FLUX:',itube,sumrec

      sumano = 0.0
      DO ic = cind1, cind2
        ip = obj(cell2obj(ic))%index(IND_FLUID)
        sumano = sumano + cell(ic)%vol * 
     .           fluid(ip,ion)%parano
      ENDDO
      IF (itube.EQ.4) 
     .  WRITE(0,*) 'FLUX:',itube,sumano



      sumflx = ABS(CalcFlux(LO,itube)) +
     .         ABS(CalcFlux(HI,itube)) 


      IF (itube.EQ.4) 
     .  WRITE(0,*) 'FLUX:',itube,sumflx

      IF (itube.EQ.4) 
     .  WRITE(0,*) 'FLUX:',itube,sumion-sumrec-sumflx


c     .                           ABS(CalcFlux(HI,itube))

c      sumflx = CalcFlux(LO,itube)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      LOGICAL FUNCTION CatchTube(a1,a2,b1,b2,itube)
      USE mod_geometry
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      INTEGER, INTENT(IN) :: itube
      REAL*8 , INTENT(IN) :: a1,a2,b1,b2

      INTEGER it,cind1,cind2,ic,iobj,isrf,ivtx(2)
      REAL*8  c1,c2,d1,d2,tab,tcd

      it = itube
      cind1 = tube(it)%cell_index(LO)
      cind2 = tube(it)%cell_index(HI)
      DO ic = cind1, cind2
        iobj = cell2obj(ic)
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
     .      tcd.GE.0.0D0.AND.tcd.LT.1.0D0) THEN
          CatchTube = .TRUE.
          EXIT
        ENDIF
      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE GenerateTargetGroups
      USE mod_geometry
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_sol28_targets
      IMPLICIT none

      INTEGER GetObject,GetTube

      INTEGER, PARAMETER :: MAXNTARGET = 100
      REAL*8 , PARAMETER :: TOL = 1.0D-06

      INTEGER itube,nlist,iobj,i,j,itarget,i1,i2,i3,i4,geometry,ixpt
      LOGICAL cont
      CHARACTER :: target_tag(4)*12
      REAL*8  x1,x2,x3,x4,y1,y2,y3,y4

      TYPE(type_target), ALLOCATABLE :: new_target(:)
      INTEGER          , ALLOCATABLE :: ilist(:),tcheck(:)

      target_tag(1) = 'inner, upper'  ! Target location specification
      target_tag(2) = 'inner, lower'
      target_tag(3) = 'outer, upper'
      target_tag(4) = 'outer, lower'

      ALLOCATE(ilist (ntube))
      ALLOCATE(tcheck(ntube))

      ALLOCATE(new_target(MAXNTARGET))

c...  Targets directly associated with the primary separatrix:      

      ntarget = 0

c...  Set the geometry:
      IF     (grid%nxpt.EQ.1.AND.grid%zxpt(1).LT.grid%z0     ) THEN
        geometry = LSND
      ELSEIF (grid%nxpt.EQ.2.AND.grid%zxpt(1).GT.grid%zxpt(2)) THEN
        geometry = UDND
      ELSE
        WRITE(0,*) 'NXPT   =',grid%nxpt
        WRITE(0,*) 'ZXPT1,2=',grid%zxpt
        WRITE(0,*) 'Z0     =',grid%z0
        STOP 'NOT TESTED YET - A'
      ENDIF

c...  A bit of a mess... but...
      DO ixpt = 1, grid%nxpt
        DO itarget = LO, HI
          nlist = 1
          IF (ixpt.EQ.1) THEN
            ilist(1) = grid%isep
          ELSE
            SELECTCASE (geometry)
              CASE (UDND)
                IF (itarget.EQ.LO) THEN
                  ilist(1) = GetTube(grid%ixpt(2,2),IND_OBJECT)
                ELSE
                  ilist(1) = GetTube(grid%ixpt(2,1),IND_OBJECT)
                ENDIF
              CASE DEFAULT
                STOP 'NOT TESTED YET - B'
            ENDSELECT
          ENDIF
          tcheck = 0
          tcheck(ilist(1)) = 1
          IF (itarget.EQ.LO) THEN
            i1 = 1
            i2 = 2
            i3 = 1
            i4 = 2
          ELSE
            i1 = 4
            i2 = 3
            i3 = 4
            i4 = 3
          ENDIF
          cont = .TRUE.
          DO WHILE (cont)
            cont = .FALSE.
            DO itube = grid%isep, ntube
              IF (tcheck(itube).NE.0) CYCLE
              iobj = GetObject(tube(itube)%cell_index(itarget),IND_CELL)  ! Replace with a function... or add to type_tube...
              CALL GetVertex(iobj,i1,x1,y1)
              CALL GetVertex(iobj,i2,x2,y2)
              DO i = nlist, 1, -1
                iobj = GetObject(tube(ilist(i))%cell_index(itarget),
     .                           IND_CELL)
                CALL GetVertex(iobj,i3,x3,y3)
                CALL GetVertex(iobj,i4,x4,y4)        
                IF     (DABS(x1-x4).LT.TOL.AND.DABS(y1-y4).LT.TOL) THEN
                  nlist = nlist + 1
                  ilist (nlist) = itube
                  tcheck(itube) = 1
                  cont = .TRUE.
                ELSEIF (DABS(x2-x3).LT.TOL.AND.DABS(y2-y3).LT.TOL) THEN
                  DO j = nlist+1, 2, -1
                    ilist(j) = ilist(j-1)
                  ENDDO
                  nlist = nlist + 1
                  ilist (1    ) = itube
                  tcheck(itube) = 1
                  cont = .TRUE.
                ENDIF
              ENDDO
            ENDDO
          ENDDO

c          WRITE(0,*) 'NLIST=',nlist
c          WRITE(0,*) 'ILIST=',ilist(1:nlist)
c          WRITE(0,*) 'TCHECK=',tcheck(1:ntube)
 
          ntarget = ntarget + 1
          IF (ntarget.GT.MAXNTARGET)
     .      CALL ER('GenerateTargetGroups','Increase MAXNTARGET',*99)

          new_target(ntarget)%version        = 1.0
          new_target(ntarget)%position       = itarget
          new_target(ntarget)%nlist          = nlist
          new_target(ntarget)%ilist(1:nlist) = ilist(1:nlist)

          SELECTCASE (geometry)
            CASE (LSND)
              IF (itarget.EQ.LO) new_target(ntarget)%location = 2
              IF (itarget.EQ.HI) new_target(ntarget)%location = 4
            CASE (USND)
              IF (itarget.EQ.LO) new_target(ntarget)%location = 3
              IF (itarget.EQ.HI) new_target(ntarget)%location = 1
            CASE (UDND)
              IF (ixpt.EQ.1) THEN
                IF (itarget.EQ.LO) new_target(ntarget)%location = 3
                IF (itarget.EQ.HI) new_target(ntarget)%location = 1
              ELSE            
                IF (itarget.EQ.LO) new_target(ntarget)%location = 2
                IF (itarget.EQ.HI) new_target(ntarget)%location = 4
              ENDIF
            CASE (LDND)
              STOP 'NOT CHECKED - A'
              IF (ixpt.EQ.1) THEN
                IF (itarget.EQ.LO) new_target(ntarget)%location = 2
                IF (itarget.EQ.HI) new_target(ntarget)%location = 4
              ELSE            
                IF (itarget.EQ.LO) new_target(ntarget)%location = 3
                IF (itarget.EQ.HI) new_target(ntarget)%location = 1
              ENDIF
            CASE (CDND)
              STOP 'NOT CHECKED - B'
              IF (ixpt.EQ.1) THEN
                IF (itarget.EQ.LO) new_target(ntarget)%location = 2
                IF (itarget.EQ.HI) new_target(ntarget)%location = 1
              ELSE            
                IF (itarget.EQ.LO) new_target(ntarget)%location = 3
                IF (itarget.EQ.HI) new_target(ntarget)%location = 4
              ENDIF
            CASE DEFAULT
              CALL ER('GenerateTargetGroups','Geometry specification '//
     .                'not recognised',*99)
          ENDSELECT
c          WRITE(0,*) 'NTARGET =',ntarget
c          WRITE(0,*) 'LOCATION=',new_target(ntarget)%location
          new_target(ntarget)%tag = 
     .      target_tag(new_target(ntarget)%location)

        ENDDO  ! target

      ENDDO  ! x-point

c...  Transfer to mod_sol28_target array:
      IF (ALLOCATED(target)) DEALLOCATE(target)
      ALLOCATE(target(ntarget))
      target(1:ntarget) = new_target(1:ntarget)
      DEALLOCATE(new_target)

c...  Add a check to make sure a tube is assigned twice, no more, no less...

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE GenerateTubeGroups
      USE mod_geometry
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      INTEGER GetObject,GetTube
      LOGICAL CatchTube
      REAL    ParticleSource

      INTEGER iobj,isrf,it,it1,it2,ic,ic1,ic2,cind1,cind2,cind3,cind4,
     .        ivtx,ivtx1,ivtx2,in,itsep,itsep2,region,i,nxpt,ixpt(2,2),
     .        region_tube(0:100,4),
     .        region_list(0:100,6)
      LOGICAL match,debug
      REAL*8  a1,a2,b1,b2,c1,c2,tab,tcd,rxpt(2),zxpt(2),r0,z0,
     .        volsrc(6)

      debug = .TRUE.

      CALL SetupCell2Obj(ncell,IND_CELL)  ! *** GET RID OF THIS! ***

c...  Calculate the magnetic axis, approximately:
      it = 1
      cind1 = tube(it)%cell_index(LO)
      cind2 = tube(it)%cell_index(HI)
      DO ic = cind1, cind2
        iobj = cell2obj(ic)
        isrf = obj(iobj)%iside(1)
        in = 1
        IF (isrf.LT.0) in = 2
        ivtx = srf(ABS(isrf))%ivtx(in)
        r0 = r0 + vtx(1,ivtx)
        z0 = z0 + vtx(2,ivtx)
      ENDDO
      r0 = r0 / DBLE(cind2 - cind1 + 1)
      z0 = z0 / DBLE(cind2 - cind1 + 1)

c...  Find x-points:
      nxpt = 0
      ixpt = 0
      rxpt = 0.0D0
      zxpt = 0.0D0

      it = grid%isep
      cind1 = tube(it)%cell_index(LO)
      cind2 = tube(it)%cell_index(HI)
      DO ic1 = cind1, cind2-2
        iobj = GetObject(ic1,IND_CELL)
        isrf = obj(iobj)%iside(1) !3)
        in = 1!2
        IF (isrf.LT.0) in = 2!1
        ivtx1 = srf(ABS(isrf))%ivtx(in)
        DO ic2 = ic1+2, cind2
          iobj = GetObject(ic2,IND_CELL)  ! cell2obj(ic2)
          isrf = obj(iobj)%iside(1)!3)
          in = 1!2
          IF (isrf.LT.0) in = 2!1
          ivtx2 = srf(ABS(isrf))%ivtx(in)
          IF (ivtx1.EQ.ivtx2) THEN  ! A bit risky this...
c            WRITE(0,*) 'ITUBE!',it
c            WRITE(0,*) 'cind1,2=',cind1,cind2
c            WRITE(0,*) 'ic1,2  =',ic1,ic2
c            WRITE(0,*) 'in     =',in
            nxpt = nxpt + 1
            rxpt(nxpt) = vtx(1,ivtx1)   
            zxpt(nxpt) = vtx(2,ivtx1)   
            ixpt(nxpt,1) = GetObject(ic1,IND_CELL)  ! cell2obj(ic1)
            ixpt(nxpt,2) = GetObject(ic2,IND_CELL)  ! cell2obj(ic2)
            itsep = it
            GOTO 10
          ENDIF
        ENDDO
      ENDDO
 10   CONTINUE

c     Check for a second x-point:
      DO it1 = grid%isep+1, ntube-1
        cind1 = tube(it1)%cell_index(LO)
        cind2 = tube(it1)%cell_index(HI)
          DO ic1 = cind1, cind2
          iobj = GetObject(ic1,IND_CELL)  ! cell2obj(ic1)
          isrf = obj(iobj)%iside(1)!3)
          in = 1!2
          IF (isrf.LT.0) in = 2!1
          ivtx1 = srf(ABS(isrf))%ivtx(in)
          DO it2 = it1+1, ntube
            cind3 = tube(it2)%cell_index(LO)
            cind4 = tube(it2)%cell_index(HI)
            DO ic2 = cind3, cind4
              iobj = GetObject(ic2,IND_CELL)  ! cell2obj(ic2)
              isrf = obj(iobj)%iside(1)!3)
              in = 1!2
              IF (isrf.LT.0) in = 2!1
              ivtx2 = srf(ABS(isrf))%ivtx(in)
              IF (ivtx1.EQ.ivtx2) THEN
c                WRITE(0,*) 'ITUBE!',it1,it2
c                WRITE(0,*) 'cind1,2=',cind1,cind2
c                WRITE(0,*) 'ic1,2  =',ic1,ic2
c                WRITE(0,*) 'in     =',in
                nxpt = nxpt + 1
                rxpt(nxpt) = vtx(1,ivtx1)   
                zxpt(nxpt) = vtx(2,ivtx1)   
                ixpt(nxpt,1) = GetObject(ic1,IND_CELL) ! cell2obj(ic1)
                ixpt(nxpt,2) = GetObject(ic2,IND_CELL) ! cell2obj(ic2)
                itsep2 = GetTube(obj(ixpt(nxpt,2))%omap(4),IND_OBJECT)
                 GOTO 20  ! Spagetti, but only a little bit
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
 20   CONTINUE

      IF (debug) THEN
        WRITE(logfp,*) 'R,Z0:',r0,z0
        WRITE(logfp,*) 'IXPT1 :',ixpt(1:nxpt,1)
        WRITE(logfp,*) 'IXPT2 :',ixpt(1:nxpt,2)
        WRITE(logfp,*) 'RXPT  :',rxpt(1:nxpt)
        WRITE(logfp,*) 'ZXPT  :',zxpt(1:nxpt)
      ENDIF

      grid%isep  = itsep
      grid%isep2 = itsep2

      grid%nxpt = nxpt
      grid%ixpt = ixpt
      grid%rxpt = SNGL(rxpt)
      grid%zxpt = SNGL(zxpt)

      region_tube = 0

c...  High field side SOL rings:
      a1 = r0
      a2 = z0
      b1 = 0.001D0
      b2 = z0
      DO it = itsep, ntube 
        IF (CatchTube(a1,a2,b1,b2,it)) THEN
          IF (debug) WRITE(logfp,*) 'TUBE CAUGHT HFS:',it
          region_tube(0,1) = region_tube(0,1) + 1
          region_tube(region_tube(0,1),1) = it
          tube(it)%type = GRD_SOL
        ENDIF
      ENDDO
c...  Identify low field side SOL rings:
      a1 = r0
      a2 = z0
      b1 = r0 + 100.0D0
      b2 = z0
      DO it = itsep, ntube 
        IF (CatchTube(a1,a2,b1,b2,it)) THEN
          IF (debug) WRITE(logfp,*) 'TUBE CAUGHT LFS:',it
          region_tube(0,2) = region_tube(0,2) + 1
          region_tube(region_tube(0,2),2) = it
          tube(it)%type = GRD_SOL
        ENDIF
      ENDDO
c...  Primary PFR:
      a1 = rxpt(1)
      a2 = zxpt(1)
      b1 = 10.0D0 * (rxpt(1) - r0) + r0
      b2 = 10.0D0 * (zxpt(1) - z0) + z0
      DO it = itsep, ntube 
        IF (CatchTube(a1,a2,b1,b2,it)) THEN
          IF (debug) WRITE(logfp,*) 'TUBE CAUGHT PRIMARY   PFR:',it
          region_tube(0,4) = region_tube(0,4) + 1
          region_tube(region_tube(0,4),4) = it
          tube(it)%type = GRD_PFZ
        ENDIF
      ENDDO

c...  Secondary PFR:
      IF (grid%nxpt.EQ.2) THEN
        a1 = rxpt(2)
        a2 = zxpt(2)
        b1 = 10.0D0 * (rxpt(2) - r0) + r0
        b2 = 10.0D0 * (zxpt(2) - z0) + z0
        DO it = itsep, ntube 
          IF (CatchTube(a1,a2,b1,b2,it)) THEN
            IF (debug) WRITE(logfp,*) 'TUBE CAUGHT SECONDARY PFR:',it
            region_tube(0,3) = region_tube(0,3) + 1
            region_tube(region_tube(0,3),3) = it
            tube(it)%type = GRD_PFZ
          ENDIF
        ENDDO
      ENDIF





      RETURN


c *** NOT EXECUTED ***

      DO it = 1, 4
        WRITE(0,*) 'REGION_TUBE:',region_tube(0:region_tube(0,it),it)
      ENDDO

c...  Associate tubes with conservation regions:
      region_list = 0
      region = 1
      DO it = 1, itsep-1
        region_list(0,region) = region_list(0,region) + 1        
        region_list(region_list(0,region),region) = it
      ENDDO
      region = 2
      DO it1 = 1, region_tube(0,1)
        match = .FALSE.
        DO it2 = 1, region_tube(0,2)        
          IF (region_tube(it1,1).EQ.region_tube(it2,2)) match = .TRUE.
        ENDDO
        IF (match) THEN
          region_list(0,region) = region_list(0,region) + 1        
          region_list(region_list(0,2),region) =  region_tube(it1,1)
        ENDIF
      ENDDO
      region = 3
      DO it1 = 1, region_tube(0,1)
        match = .FALSE.
        DO it2 = 1, region_tube(0,2)        
          IF (region_tube(it1,1).EQ.region_tube(it2,2)) match = .TRUE.
        ENDDO
        IF (.NOT.match) THEN
          region_list(0,region) = region_list(0,region) + 1        
          region_list(region_list(0,region),region) = region_tube(it1,1)
        ENDIF
      ENDDO
      region = 4
      DO it1 = 1, region_tube(0,2)
        match = .FALSE.
        DO it2 = 1, region_tube(0,1)        
          IF (region_tube(it1,2).EQ.region_tube(it2,1)) match = .TRUE.
        ENDDO
        IF (.NOT.match) THEN
          region_list(0,region) = region_list(0,region) + 1        
          region_list(region_list(0,region),region) = region_tube(it1,2)
        ENDIF
      ENDDO

      region_list(:,5) = region_tube(:,3)
      region_list(:,6) = region_tube(:,4)

c...  Add up total ionisation source in the core (need to subtract rate of increase
c     of core inventory and the particles pumped by the core inner boundary...):
      volsrc = 0.0D0
      region = 1
      DO i = 1, region_list(0,region)
        it = region_list(i,region)

        it = 4
        volsrc(region) = volsrc(region) + ParticleSource(it)
              STOP 'sdgsd --- b'
      ENDDO


      DO it = 1, 6
        WRITE(0,*) 'REGION_LIST:',region_list(0:region_list(0,it),it)
      ENDDO





!...  Need to put this in mod_geometry.f:
      DEALLOCATE(cell2obj)

      STOP 'sdgsdgs'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadMaterialData
      USE mod_sol28_global
      IMPLICIT none

      INTEGER   matfp
      CHARACTER buffer*1024

      matfp = 99

      nmaterial = 0

      SELECTCASE(opt%mat_opt) 
        CASE(1)  
          OPEN(UNIT=matfp,FILE=TRIM(opt%mat_file),ACCESS='SEQUENTIAL',
     .         ERR=98)     
          DO WHILE(.TRUE.) 
            READ(matfp,'(A1024)',END=10) buffer        
            IF (buffer(1:1).EQ.'*'.OR.LEN_TRIM(buffer).LT.10) CYCLE  ! A bit arbitrary on the length test...
            nmaterial = nmaterial + 1
            READ(buffer,*) 
     .        material(nmaterial)%tag,
     .        material(nmaterial)%type,
     .        material(nmaterial)%A,
     .        material(nmaterial)%Z,
     .        material(nmaterial)%mass
          ENDDO
 10       CONTINUE
          CLOSE(matfp)
        CASE DEFAULT
          CALL ER('LoadMaterialData','Unknown format',*99)
      ENDSELECT

      RETURN
 98   CALL ER('LoadMaterialData','File access error',*99)
      WRITE(0,*) '  FILE= ',TRIM(opt%mat_file)
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE SequenceWall(mode)
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_sol28_wall
      IMPLICIT none

      INTEGER, INTENT(IN) :: mode

      REAL*8 , PARAMETER :: DTOL = 1.0D-06

      INTEGER i1,i2,nwall_1
      LOGICAL cont
      TYPE(type_wall) :: tmp_wall


c...  Put all of the CLASS=1 (main continuous wall) to the beginning
c     of the wall segment list:
      DO i1 = 1, nwall-1
        IF (wall(i1)%class.EQ.1) CYCLE
        DO i2 = i1+1, nwall
          IF (wall(i2)%class.EQ.1) THEN
            wall(nwall+1) = wall(i1     )
            wall(i1     ) = wall(i2     )
            wall(i2     ) = wall(nwall+1)
          ENDIF
        ENDDO
      ENDDO

      DO i1 = 1, nwall
        IF (wall(i1)%class.EQ.1) nwall_1 = i1
      ENDDO
c      WRITE(0,*) 'nwall,_1:',nwall,nwall_1

c...  Re-index the wall so that it starts at the outer edge of
c     the first low index target segment, as in DIVIMP:
      SELECTCASE (mode)
        CASE (1)
        CASE (2)
          DO i1 = 1, nwall_1
            IF (wall(i1)%index(WAL_TUBE  ).EQ.grid%isep.AND.
     .          wall(i1)%index(WAL_TARGET).EQ.LO) EXIT
          ENDDO
c          WRITE(0,*) 'SEP:',i1
          IF (i1.EQ.nwall+1) 
     .      CALL ER('SequenceWall','Separatrix not found',*99)
          cont = .TRUE.
          DO WHILE (cont)
            cont = .FALSE.
            DO i2 = 1, nwall_1
              IF (DABS(wall(i2)%v1(1)-wall(i1)%v2(1)).LT.DTOL.AND.
     .            DABS(wall(i2)%v1(2)-wall(i1)%v2(2)).LT.DTOL) THEN
                IF (wall(i2)%index(WAL_TARGET).EQ.LO) cont = .TRUE.
                i1 = i2
                EXIT
              ENDIF
            ENDDO
            IF (i2.EQ.nwall_1+1) 
     .        CALL ER('SequenceWall','Origin not set',*99)
          ENDDO
c         Swap the first segment with the idenfied one:           
          tmp_wall = wall(1)
          wall(1)  = wall(i1)
          wall(i1) = tmp_wall
        CASE DEFAULT
         CALL ER('SequenceWall','Unknown MODE',*99)
      ENDSELECT

c...  Clockwise order the segments that make up the standard
c     continuous wall:
      DO i1 = 1, nwall_1-1
        DO i2 = i1+1, nwall_1
          IF (DABS(wall(i1)%v2(1)-wall(i2)%v1(1)).LT.DTOL.AND.
     .        DABS(wall(i1)%v2(2)-wall(i2)%v1(2)).LT.DTOL) THEN
            IF (i2.GT.i1+1) THEN
              wall(nwall+1) = wall(i1+1   )
              wall(i1+1   ) = wall(i2     )
              wall(i2     ) = wall(nwall+1)
            ENDIF
c...        Make sure the wall closes exactly:
            IF (wall(i1  )%index(WAL_TARGET).EQ.0.AND.
     .          wall(i1+1)%index(WAL_TARGET).NE.0) THEN
              wall(i1  )%v2 = wall(i1+1)%v1
            ELSE
              wall(i1+1)%v1 = wall(i1  )%v2
            ENDIF
            EXIT
          ENDIF
        ENDDO
        IF (i2.EQ.nwall_1+1) 
     .    CALL ER('SequenceWall','The standard wall is not '//
     .            'continuous',*99)
      ENDDO
c...  Final test:
      IF (DABS(wall(nwall_1)%v2(1)-wall(1)%v1(1)).GT.DTOL.OR.
     .    DABS(wall(nwall_1)%v2(2)-wall(1)%v1(2)).GT.DTOL) 
     .  CALL ER('SequenceWall','The standard wall does not close '//
     .          'on itself',*99)
      wall(nwall_1)%v2 = wall(1)%v1

      RETURN
 99   WRITE(0,*) ' MODE = ',mode
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE MapTargetToWall(index,itube,itarget,tmp_nwall,tmp_wall)
      USE mod_geometry
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_sol28_wall
      IMPLICIT none

      INTEGER, INTENT(IN)  :: itube,itarget,tmp_nwall
      INTEGER, INTENT(OUT) :: index
      TYPE(type_wall) :: tmp_wall(tmp_nwall)

      INTEGER GetObject

      REAL*8 , PARAMETER :: DTOL=1.0D-06

      INTEGER iwall,fp,iobj
      LOGICAL debug
      REAL*8  x1,x2,x3,x4,y1,y2,y3,y4,s12,s34,s12max,length

      debug = .TRUE.
      fp = 88

      index = -1

      iobj = GetObject(tube(itube)%cell_index(itarget),IND_CELL)                   

      IF (debug) WRITE(fp,*) 'WALL-TARGET INTERSECTION',
     .                       itube,iobj,itarget

      CALL GetVertex(iobj,3,x3,y3)
      CALL GetVertex(iobj,4,x4,y4)
      x1 = 0.5D0 * (x3 + x4)
      y1 = 0.5D0 * (y3 + y4)
      CALL GetVertex(iobj,1,x3,y3)
      CALL GetVertex(iobj,2,x4,y4)
      x2 = 0.5D0 * (x3 + x4)
      y2 = 0.5D0 * (y3 + y4)
      IF (itarget.EQ.HI) THEN
        x3 = x1
        y3 = y1
        x1 = x2
        y1 = y2
        x2 = x3
        y2 = y3
      ENDIF
c     Increase the length of the line segment forward and backward 
c     in case there's a nominal mismatch between the wall and target 
c     specifications or if the line segment is very short:
      length = DSQRT((x1-x2)**2 + (y1-y2)**2)
      x2 = x1 + MAX(2.0D0,0.1D0 / length) * (x2 - x1)
      y2 = y1 + MAX(2.0D0,0.1D0 / length) * (y2 - y1)
      length = DSQRT((x1-x2)**2 + (y1-y2)**2)
      x1 = x2 + MAX(2.0D0, 0.2D0 / length) * (x1 - x2)
      y1 = y2 + MAX(2.0D0, 0.2D0 / length) * (y1 - y2)

c     Search the wall for intersections:
      s12max = 1.0D+10
      DO iwall = 1, tmp_nwall
        x3 = tmp_wall(iwall)%v1(1)
        y3 = tmp_wall(iwall)%v1(2)
        x4 = tmp_wall(iwall)%v2(1)
        y4 = tmp_wall(iwall)%v2(2)
        CALL CalcInter(x1,y1,x2,y2,x3,y3,x4,y4,s12,s34) 
        IF (debug) THEN
          WRITE(fp,*) ' --------------------'
          WRITE(fp,*) ' X,Y1=',x1,y1
          WRITE(fp,*) ' X,Y2=',x2,y2
          WRITE(fp,*) '  IWALL    :',iwall
          WRITE(fp,*) '  S12,S34  :',s12,s34
          WRITE(fp,*) '  X,Y3     :',x3,y3
          WRITE(fp,*) '  X,Y4     :',x4,y4
          IF (s34.GT.0.0D0.AND.s34.LE.1.0D0) 
     .      WRITE(fp,*) '     **** BING **** '
        ENDIF
        IF (s12.GT.0.0D0.AND.s12.LT.1.0D0.AND.
     .      s34.GT.0.0D0.AND.s34.LE.1.0D0.AND.
     .      s12.LT.s12max) THEN
          s12max = s12
          index = iwall
          IF (debug) WRITE(fp,*) '  *** INTERSECTION ***',s12,iwall
        ENDIF
      ENDDO


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE ProcessWall
      USE mod_geometry
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_sol28_wall
      IMPLICIT none


      INTEGER, PARAMETER :: MAXNTMPWALL = 10000

      INTEGER GetObject

      INTEGER   fp,walfp,iopt,i1,i2,i3,itube,iobj,itarget,nwall_1,iwall,
     .          tmp_nwall,imat
      CHARACTER buffer*1024

      REAL*8         , ALLOCATABLE :: rwall(:,:),zwall(:,:)
      TYPE(type_wall), ALLOCATABLE :: tmp_wall(:)

      IF (ALLOCATED(wall)) DEALLOCATE(wall)
      ALLOCATE(wall(MAXNTMPWALL))

c...  Build the list of line segments that make up the wall:   
      nwall = 0
      DO iopt = 1, nopt_wall
c        WRITE(0,*) 'NOPT_WALL:',iopt,nopt_wall
        SELECTCASE(opt_wall(iopt)%file_format)
          CASE(1)
            walfp = 99
            OPEN(UNIT=walfp,FILE=TRIM(opt_wall(iopt)%file_name),
     .           ACCESS='SEQUENTIAL',ERR=98)     
            DO WHILE(.TRUE.) 
              READ(walfp,'(A1024)',END=10) buffer        
              IF (buffer(1:1).EQ.'*'.OR.LEN_TRIM(buffer).LT.5) CYCLE  ! A bit arbitrary on the length test...
              nwall = nwall + 1
              wall(nwall) = opt_wall(iopt)
              READ(buffer,*) wall(nwall)%v1(1:2),wall(nwall)%v2(1:2)
c...          Defaults:
              wall(nwall)%index(WAL_INDEX ) = 0
              wall(nwall)%index(WAL_TUBE  ) = 0
              wall(nwall)%index(WAL_TARGET) = 0
              wall(nwall)%material          = 0
c...          Check if end points are degenerate:
              IF (wall(nwall)%v1(1).EQ.wall(nwall)%v2(1).AND.
     .            wall(nwall)%v1(2).EQ.wall(nwall)%v2(2)) THEN
                CALL WN('ProcessWall','Degenerate end points '//
     .                  'identified and discarded',*99)
                WRITE(0,*) ' FNAME= ',TRIM(opt_wall(iopt)%file_name)
                WRITE(0,*) ' NWALL= ',nwall
                nwall = nwall - 1
              ENDIF
            ENDDO
 10         CONTINUE
            CLOSE(walfp)
          CASE DEFAULT
            CALL ER('ProcessWall','Unknown format',*99)
        ENDSELECT
      ENDDO
c      WRITE(0,*) 'NWALL:',nwall

      CALL SaveWallGeometry

c...  Sequence the wall so that each line segment in the list follows 
c     the once it's geometrically connected to:
c      WRITE(0,*) 'SEQUENCE WALL ONCE'
      CALL SequenceWall(1)
      

c      CALL GenerateOutputFiles
c      STOP 'sdfsdfsd'

c...  Store the pre-clipped wall so that it can be used to determine
c     the properties of the target segments:
      ALLOCATE(tmp_wall(nwall))
      tmp_nwall = nwall
      tmp_wall  = wall

c...  Clip the standard wall to the grid:
      DO i1 = 1, nwall
        IF (wall(i1)%class.EQ.1) nwall_1 = i1
      ENDDO
      ALLOCATE(rwall(nwall_1,2))
      ALLOCATE(zwall(nwall_1,2))
      rwall(1:nwall_1,1) = wall(1:nwall_1)%v1(1)  
      zwall(1:nwall_1,1) = wall(1:nwall_1)%v1(2)  
      rwall(1:nwall_1,2) = wall(1:nwall_1)%v2(1)  
      zwall(1:nwall_1,2) = wall(1:nwall_1)%v2(2)  
      CALL ClipWallToGrid(nwall_1,rwall,zwall,NWALL_1,.FALSE.)
      wall(1:nwall_1)%v1(1) = rwall(1:nwall_1,1)
      wall(1:nwall_1)%v1(2) = zwall(1:nwall_1,1)
      wall(1:nwall_1)%v2(1) = rwall(1:nwall_1,2)
      wall(1:nwall_1)%v2(2) = zwall(1:nwall_1,2)
      DEALLOCATE(rwall)
      DEALLOCATE(zwall)

c...  Delete marked segments:
      DO i1 = nwall, 1, -1
        IF (wall(i1)%v1(1).EQ.-9.99D0) THEN
          DO i2 = i1, nwall-1
            wall(i2) = wall(i2+1)
          ENDDO
          nwall = nwall - 1
        ENDIF
      ENDDO

c...  Add target segments:
      DO itube = 1, ntube
        IF (itube.GE.grid%isep) THEN
          DO itarget = LO, HI
            IF (itarget.EQ.LO) THEN
              i1 = 1
              i2 = 2
            ELSE
              i1 = 3
              i2 = 4
            ENDIF
            iobj = GetObject(tube(itube)%cell_index(itarget),IND_CELL)                   
            nwall = nwall + 1
            wall(nwall)%index = 0
            wall(nwall)%type              = -1
            wall(nwall)%class             =  1
            wall(nwall)%material_tag      = '?'
            wall(nwall)%material          = 0
            wall(nwall)%temperature       = 300.0
            wall(nwall)%index(WAL_GROUP)  = -1
            wall(nwall)%index(WAL_TUBE  ) = itube
            wall(nwall)%index(WAL_TARGET) = itarget
            CALL GetVertex(iobj,i1,wall(nwall)%v1(1),wall(nwall)%v1(2))
            CALL GetVertex(iobj,i2,wall(nwall)%v2(1),wall(nwall)%v2(2))
          ENDDO
        ENDIF
      ENDDO

c...  Check TMP_WALL to see what the surface properties of the target
c     segments should be:
      DO iwall = 1, nwall
        IF (wall(iwall)%index(WAL_TARGET).EQ.0) CYCLE
        itube   = wall(iwall)%index(WAL_TUBE  )
        itarget = wall(iwall)%index(WAL_TARGET)
        CALL MapTargetToWall(i1,itube,itarget,tmp_nwall,tmp_wall)
        IF (i1.EQ.-1) THEN
          CALL ER('ProcessWall','Target mapping failed',*99)
        ELSE
          wall(iwall)%index(WAL_GROUP) = tmp_wall(i1)%index(WAL_GROUP)
          wall(iwall)%material_tag = tmp_wall(i1)%material_tag
          wall(iwall)%temperature  = tmp_wall(i1)%temperature
        ENDIF
      ENDDO
      DEALLOCATE(tmp_wall)   

c...  Sequence the standard wall once again:
c      WRITE(0,*) 'SEQUENCE WALL TWICE'
      CALL SequenceWall(2)

c...  Resize the wall array so memory isn't wasted (assuming this procedure
c     doesn't fragment the memory space, wasting space...):
      ALLOCATE(tmp_wall(nwall))
      tmp_wall = wall
      DEALLOCATE(wall)
      ALLOCATE(wall(nwall))
      wall = tmp_wall
      DEALLOCATE(tmp_wall)

c...  Assign indexing within each class of wall segment (CLASS=1 is the
c     standard continuous wall, CLASS=2 are additional wall surfaces,
c     CLASS=3 are holes in the EIRENE triangle mesh):
      DO i1 = 1, 100  ! Assume there are never more than 100 classes...
        i2 = 0
        DO iwall = 1, nwall
          IF (wall(iwall)%class.EQ.i1) THEN
            i2 = i2 + 1
            wall(iwall)%index(WAL_INDEX) = i2
          ENDIF
        ENDDO
      ENDDO

c...  Assign wall to the appropriate entry in the list of defined materials:
      DO iwall = 1, nwall      
        DO imat = 1, nmaterial
          IF (TRIM(wall    (iwall)%material_tag).EQ.
     .        TRIM(material(imat )%tag)) THEN
            wall(iwall)%material = imat
            EXIT
          ENDIF
        ENDDO
      ENDDO


      CALL SaveWallGeometry

      CALL DumpData_OSM('output.wall_load','Done loading the wall')      

      RETURN
 98   CALL ER('ProcessWall','File access error',*99)
 99   WRITE(0,*) 'IOPT=  ',iopt
      WRITE(0,*) 'FILE=  ',TRIM(opt_wall(iopt)%file_name)
      WRITE(0,*) 'FORMAT=',opt_wall(iopt)%file_format
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadSupplimentalGridData
      USE mod_geometry 
      USE mod_grid
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      INTEGER GetObject

      INTEGER   fp,ipos,itube,version,ierr,icell,iobj,maxitube,
     .          maxipos(ntube),hold_itube
      LOGICAL   first_try,debug
      CHARACTER buffer*1024
      REAL      rdum1,rdum2
      REAL*8    a1,a2,b1,b2

      debug = .FALSE.

      IF (log.GT.0) THEN
        WRITE(logfp,*)
        WRITE(logfp,*) 'LOADING SUPPLIMENTAL GRID DATA'
      ENDIF

      first_try = .TRUE.
      fp = 99

c...  Look for <gridname>.sup:
 10   OPEN(UNIT=fp,FILE=TRIM(opt%f_grid_file)//'.sup',
     .     ACCESS='SEQUENTIAL',STATUS='OLD',IOSTAT=ierr)     
      IF (debug) WRITE(0,*) 'IERR=',ierr

c...  Check if the supplimental data file is not around, and if not,
c     call the requisite IDL routines to generate the data:
      IF (ierr.NE.0) THEN
        IF (.NOT.first_try) 
     .    CALL ER('ProcessGrid','Trouble generating '//
     .            'the supplimental grid data',*99)
c       Output the OSM geometry data file that's read in by IDL:
        OPEN (UNIT=fp,FILE='grid.sup',ACCESS='SEQUENTIAL',STATUS='NEW')      
        WRITE(fp,'(A)') '* Grid geometry file from OSM'
        WRITE(fp,'(A,1X,2A6,6A12)') 
     .    '*','IP','IT','Rcen (m)','Zcen (m)','Rmid1 (m)',
     .    'Zmid1 (m)','R1 (m)','Z1 (m)'
        WRITE(fp,*) 3  ! Format code
        DO itube = 1, ntube
          DO icell = tube(itube)%cell_index(LO),
     .               tube(itube)%cell_index(HI)
            iobj = GetObject(icell,IND_CELL)                   
            CALL GetVertex(iobj,1,a1,b1)
            CALL GetVertex(iobj,2,a2,b2)
            WRITE(fp,'(2X,2I6,6F12.7)') 
     .        icell-tube(itube)%cell_index(LO)+1,itube,
     .        cell(icell)%cencar(1),cell(icell)%cencar(2),
     .        0.5*(a1+a2),0.5*(b1+b2),a1,b1
          ENDDO
        ENDDO
        CLOSE(fp)
c       Call IDL:
        CALL CIssue('idl grid_run.pro -quiet -args suppliment '//
     .              'grid.sup '//TRIM(opt%f_grid_file)//'.equ',ierr)
c       Copy the data to <gridname>.sup for storage in the equilibrium 
c       directory (the file is moved by the OSM run script):
        CALL CIssue('cp grid.sup.out '//TRIM(opt%f_grid_file)//
     .              '.sup',ierr)
c       Now try to open the supplimental data file:
        first_try = .FALSE.
        GOTO 10  ! Spagetti...
      ENDIF

c...  Load the supplimental data from the IDL file:
      tube(1:ntube)%psin = 0.0
      hold_itube = -1
      maxitube = 0
      maxipos  = 0
      READ(fp,*) version  ! First line 
      DO WHILE (.TRUE.) 
        READ(fp,'(A)',END=20) buffer 
        IF (buffer(1:1).EQ.'*'.OR.buffer(1:1).EQ.'$') CYCLE  ! Comment line indicators
        SELECTCASE (version)
          CASE (1)
            READ(buffer,*,ERR=98) ipos,itube,rdum1,rdum2
            maxitube = MAX(maxitube,itube)
            IF (maxitube.GT.ntube) EXIT
            maxipos(itube) = ipos
            tube(itube)%psin = tube(itube)%psin + rdum1  ! Take an average PSIn value for the tube
            IF (itube.NE.hold_itube) THEN
              hold_itube = itube
              icell = tube(itube)%cell_index(LO) - 1
            ENDIF
            icell = icell + 1
            field(icell)%bratio = rdum2
          CASE DEFAULT
            CALL ER('ProcessGrid','Unrecognised version number for '//
     .              'the supplimental grid data file',*99)
        ENDSELECT
      ENDDO
 20   CLOSE(fp)
      DO itube = 1, ntube
        IF (maxipos(itube).NE.tube(itube)%n.OR.maxitube.NE.ntube) THEN 
c         There may have been some modifications to the grid which
c         are case dependent, so rebuild the supplimental file.  This 
c         isn't the best solution perhaps but it avoids having to save
c         endless sub-versions of the file, one for each case:
          IF (first_try) THEN
            CALL WN('ProcessGrid','Supplimental data file not '//
     .              'consistent with the size of the grid, '// 
     .              'rebuilding the file')
            CALL CIssue('rm -f grid.sup',ierr)
            CALL CIssue('rm -f '//TRIM(opt%f_grid_file)//'.sup',ierr)
            GOTO 10  ! Sorry for the spagetti...
          ELSE
            CALL ER('ProcessGrid','Supplimental data file not '//
     .              'consistent with the size of the grid',*99) 
          ENDIF
        ELSE
          tube(itube)%psin = tube(itube)%psin / REAL(tube(itube)%n)
        ENDIF
      ENDDO

      IF (.TRUE.) THEN
        DO itube = 1, ntube
         tube(itube)%bratio(LO)=field(tube(itube)%cell_index(LO))%bratio
         tube(itube)%bratio(HI)=field(tube(itube)%cell_index(HI))%bratio
        ENDDO
        CALL osm_DeriveGridQuantities
      ENDIF

      IF (log.GT.0) WRITE(logfp,*) 'DONE'

      RETURN
 98   CALL ER('ProcessGrid','Format error in supplimental file',*99)
 99   WRITE(0,*) 'IPOS,ITUBE  =',ipos,itube
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE osm_DeriveGridQuantities
      USE mod_geometry 
      USE mod_grid
      USE mod_legacy
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_solps
      IMPLICIT none

      INTEGER GetObject
      REAL*8  CalcPolygonArea

      INTEGER itube,ic1,ic2,icell,iobj,i
      REAL    deltap,deltas,area,volume
      REAL*8  a(3),b(3),c(3),d(3),e(3),f(3),x(10),y(10)

      DO itube = 1, ntube
        ic1 = tube(itube)%cell_index(LO)
        ic2 = tube(itube)%cell_index(HI)
        cell(ic1)%sbnd = 0.0
        cell(ic1)%pbnd = 0.0

        DO icell = ic1, ic2
          iobj = GetObject(icell,IND_CELL)

          CALL GetVertex(iobj,1,a(1),a(2))
          CALL GetVertex(iobj,2,b(1),b(2))
          CALL GetVertex(iobj,3,c(1),c(2))
          CALL GetVertex(iobj,4,d(1),d(2))

          e(:) = 0.5D0 * (a(:) + b(:))
          f(:) = 0.5D0 * (c(:) + d(:))

          deltap = SNGL(DSQRT((e(1)-f(1))**2 + (e(2)-f(2))**2))
          deltas = deltap / (field(icell)%bratio + EPS10)

          IF (icell.GT.ic1) THEN
            cell(icell)%pbnd(1) = cell(icell-1)%pbnd(2)
            cell(icell)%sbnd(1) = cell(icell-1)%sbnd(2)
          ENDIF
          cell(icell)%p  = cell(icell)%pbnd(1) + 0.5 * deltap
          cell(icell)%s  = cell(icell)%sbnd(1) + 0.5 * deltas
          cell(icell)%ds = deltas
          cell(icell)%pbnd(2) = cell(icell)%pbnd(1) + deltap
          cell(icell)%sbnd(2) = cell(icell)%sbnd(1) + deltas

c...      Calculate cell area and volume:
          DO i = 1, obj(iobj)%nside
            CALL GetVertex(iobj,i,x(i),y(i))
          ENDDO
          area   = SNGL(CalcPolygonArea(x,y,obj(iobj)%nside))
          volume = 2.0 * V_PI * cell(icell)%cencar(1) * area
          cell(icell)%vol = volume
        ENDDO

        tube(itube)%pmax = cell(ic2)%pbnd(2)
        tube(itube)%smax = cell(ic2)%sbnd(2)
      ENDDO
 
      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE ProcessGrid
      USE mod_geometry 
      USE mod_grid
c      USE mod_grid_divimp
      USE mod_legacy
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_solps
      IMPLICIT none

      INTEGER CalcPoint
      REAL*8  CalcPolygonArea

      INTEGER   status, grd_load_method,ir,ik,iside,i1,i2,id,i,target,
     .          fp,itube,icell,iobj
      LOGICAL   debug
      REAL      deltap,deltas,area,volume
      REAL*8    a(3),b(3),c(2,3),d(2,3),e(3),f(3),g(2,3),t(2),
     .          a1,a2,b1,b2,c1,c2,shyp,sadj,x(10),y(10)
      TYPE(type_srf   ) newsrf
      TYPE(type_object) newobj

      INTEGER, PARAMETER :: GRD_LOAD_NEW = 2, GRD_LOAD_OLD = 1

      debug = .FALSE.

      grd_load_method = 2

      SELECTCASE (grd_load_method)
        CASE (GRD_LOAD_NEW)
          CALL LoadGeneralisedGrid
c...      Assign geometry structures:
          irsep      = grid_load%irsep      
          irsep2     = grid_load%irsep2     
          irwall     = grid_load%irwall     
          irtrap     = grid_load%irtrap     
          nrs        = grid_load%nrs        
          ikti       = grid_load%ikti       
          ikto       = grid_load%ikto       
          nks(1:nrs) = grid_load%nks(1:nrs)    
          IF (debug) THEN
            WRITE(0,*) nrs        
            WRITE(0,*) nknot
            WRITE(0,*) irsep      
            WRITE(0,*) irsep2     
            WRITE(0,*) irwall     
            WRITE(0,*) irtrap     
            WRITE(0,*) ikti       
            WRITE(0,*) ikto       
            WRITE(0,*) nks(1:nrs)
          ENDIF
c...      Define arrays:
          ntube     = nrs
          ncell     = nknot
          nion      = 1
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
          ALLOCATE(pin     (npin     ,nion))   ! *** REDUNDANT -- SHOULD USE kinetic? ***  Maybe not, since arrays may be different, ie one 3D, the other 2D?
          ALLOCATE(photon  (nphoton  ,nion))   
          ALLOCATE(drift   (ndrift   ,nion))
          ALLOCATE(kinetic (nkinetic ,nion))  
          ALLOCATE(fluid   (nfluid   ,nion))
          ALLOCATE(impurity(nimpurity,nion))   ! *** REDUNDANT -- SHOULD USE kinetic? ***
c...      Group association:
          ngrp = 1
          grp(1)%origin = GRP_MAGNETIC_GRID  ! *** NEED AddGroup and CopyGroup functions ***
          grp(1)%type   = GRP_QUADRANGLE 
c...      Set grid quantities:
          grid%n     = nrs
          grid%isep  = irsep   ! These are over-written in GenerateTubeGroups
          grid%isep2 = irsep2  !   (this one)
          grid%ipfz  = irtrap 
          grid%ikto  = ikto
          grid%ikti  = ikti
c...      Setup for SOLPS data mapping:
c          IF (solps_opt.GT.0) THEN
c            divimp_maxnrs = nrs + 3
c            divimp_maxnks = 0
c            DO ir = 1, nrs
c              divimp_maxnks = MAX(divimp_maxnks,nks(ir))
c            ENDDO
c            CALL ALLOC_GRID(divimp_maxnks,divimp_maxnrs)
c          ENDIF
c...      Assign TUBE and CELL quanitites:
          nobj = 0
          nsrf = 0
          nvtx = 0
          ncell = 0
          DO ir = 1, nrs
            tube(ir)%bratio = 0.0
            tube(ir)%dds    = 0.0
            tube(ir)%rp     = 0.0
            tube(ir)%costet = 0.0
c...
            tube(ir)%type = GRD_PFZ
            IF (ir.LE.irwall) tube(ir)%type = GRD_SOL
            IF (ir.LT.irsep ) tube(ir)%type = GRD_CORE
            tube(ir)%ir = ir + 1
            IF (ir.GT.irtrap) tube(ir)%ir = tube(ir)%ir + 2
            tube(ir)%n = nks(ir)
            IF (ir.EQ.1) THEN
              tube(ir)%cell_index(1) = 1
            ELSE
              tube(ir)%cell_index(1) = tube(ir-1)%cell_index(2) + 1
            ENDIF
            tube(ir)%cell_index(2) = tube(ir)%cell_index(1) + nks(ir)-1
c...
            cell(ncell+1)%pbnd = 0.0
            cell(ncell+1)%sbnd = 0.0
            DO ik = 1, nks(ir)
              id = imap(ik,ir)
              ncell = ncell + 1
              cell(ncell)%ik        = ik
              cell(ncell)%ir        = tube(ir)%ir
              cell(ncell)%nside     = 4
              cell(ncell)%cencar(1) = knot(id)%rcen
              cell(ncell)%cencar(2) = knot(id)%zcen
              cell(ncell)%cencar(3) = 0.0D0
c...          For SOLPS data mapping:
c              IF (ALLOCATED(divimp_ik)) THEN
c                divimp_ik(ik,ir) = knot(id)%ik 
c                divimp_ir(ik,ir) = knot(id)%ir
c              ENDIF

              field(ncell)%bratio = knot(id)%bratio

              newobj%group         = ngrp
              newobj%index(IND_IK) = ik
              newobj%index(IND_IR) = ir
              newobj%index(IND_IS) = 0
              newobj%index(IND_CELL   ) = ncell
              newobj%index(IND_FLUID  ) = ncell
              newobj%index(IND_KINETIC) = ncell
              newobj%index(IND_NEUTRAL) = ncell
              newobj%index(IND_FIELD  ) = ncell
              newobj%segment(1) = 0
              newobj%phi        = 0.0
              newobj%nside      = 4
              DO iside = 1, 4
                newsrf%type = SPR_LINE_SEGMENT
                newsrf%obj  = ncell
                newsrf%side = iside
                newsrf%nvtx = 2
                i1 = iside
                i2 = iside + 1
                IF (i2.EQ.5) i2 = 1
                a(1) = knot(id)%rv(i1)
                a(2) = knot(id)%zv(i1)
                a(3) = 0.0D0
                newsrf%ivtx(1) = AddVertex(a) 
                b(1) = knot(id)%rv(i2)
                b(2) = knot(id)%zv(i2)
                b(3) = 0.0D0
                newsrf%ivtx(2) = AddVertex(b) 
                newobj%iside(iside) = AddSurface(newsrf)
                IF (iside.EQ.1) THEN
                  c(1,:) = a(:)
                  c(2,:) = b(:)
                ENDIF
                IF (iside.EQ.3) THEN
                  d(1,:) = a(:)
                  d(2,:) = b(:)
                ENDIF
              ENDDO
              nobj = AddObject(newobj)
          
              e(:) = 0.5D0 * (c(1,:) + c(2,:))
              f(:) = 0.5D0 * (d(1,:) + d(2,:))

c Moved to osm_DeriveGridQuantities...
c              deltap = SNGL(DSQRT((e(1)-f(1))**2 + (e(2)-f(2))**2))
c              deltas = deltap / (field(ncell)%bratio + EPS10)
c              IF (ik.GT.1) THEN
c                cell(ncell)%pbnd(1) = cell(ncell-1)%pbnd(2)
c                cell(ncell)%sbnd(1) = cell(ncell-1)%sbnd(2)
c              ENDIF
c              cell(ncell)%p  = cell(ncell)%pbnd(1) + 0.5 * deltap
c              cell(ncell)%s  = cell(ncell)%sbnd(1) + 0.5 * deltas
c              cell(ncell)%ds = deltas
c              cell(ncell)%pbnd(2) = cell(ncell)%pbnd(1) + deltap
c              cell(ncell)%sbnd(2) = cell(ncell)%sbnd(1) + deltas
cc...          Calculate cell area and volume:
c              DO i = 1, obj(nobj)%nside
c                CALL GetVertex(nobj,i,x(i),y(i))
c              ENDDO
c              area   = SNGL(CalcPolygonArea(x,y,obj(nobj)%nside))
c              volume = 2.0 * V_PI * cell(ncell)%cencar(1) * area
c              cell(ncell)%vol = SNGL(volume)

c...          Calculate target quantities:
              target = 0
              IF (ir.GE.irsep) THEN
                IF     (ik.EQ.1      ) THEN
                  target = LO
                  a1 = e(1)
                  a2 = e(2)
                  g = c
                ELSEIF (ik.EQ.nks(ir)) THEN
                  target = HI
                  a1 = f(1)
                  a2 = f(2)
                  g = d
                ENDIF
              ENDIF
              IF (target.GT.0) THEN
                b1 = knot(id)%rcen
                b2 = knot(id)%zcen
                DO i = 1, 2
                  c1 = g(i,1)
                  c2 = g(i,2)
                  IF (CalcPoint(a1,a2,b1,b2,c1,c2,t(i)).LT.0) 
     .              CALL ER('ProcessGrid','Intersection not found',*99)
                ENDDO
                IF (t(1).GT.t(2)) THEN
                  c1 = g(1,1)
                  c2 = g(1,2)
                ELSE
                  t(1) = t(2)
                  c1 = g(2,1)
                  c2 = g(2,2)
                ENDIF
                b1 = a1 + t(1) * (b1 - a1)
                b2 = a2 + t(1) * (b2 - a2)
                shyp = DSQRT((a1 - c1)**2 + (a2 - c2)**2)
                sadj = DSQRT((b1 - c1)**2 + (b2 - c2)**2)

                tube(ir)%bratio(target)=field(ncell)%bratio
                tube(ir)%dds   (target)=SNGL(DSQRT((g(1,1)-g(2,1))**2 +
     .                                             (g(1,2)-g(2,2))**2))
                tube(ir)%rp    (target)=SNGL(a1)
                tube(ir)%costet(target)=SNGL(sadj / (shyp + DPS10))
c                IF (ir.EQ.13) THEN
c                  WRITE(0,*) '--------------------'
c                  WRITE(0,*) 'IR     =',ir
c                  WRITE(0,*) 'target =',target
c                  WRITE(0,*) 'C      =',c(1,1:2)
c                  WRITE(0,*) 'D      =',d(1,1:2)
c                  WRITE(0,*) 'G1     =',g(1,1:2)
c                  WRITE(0,*) 'G2     =',g(2,1:2)
c                  WRITE(0,*) 'DDS    =', tube(ir)%dds(target)
c                ENDIF
              ENDIF
            ENDDO ! End of IK loop

            tube(ir)%pmax = cell(ncell)%pbnd(2)
            tube(ir)%smax = cell(ncell)%sbnd(2)
          ENDDO
c...      Store global geometry data:   *** NOT ACCURATE, GET FROM EQUILIBRIUM! ***
          i1 = tube(1)%cell_index(1)
          i2 = tube(1)%cell_index(2)
          grid%r0 = DBLE(SUM(cell(i1:i2)%cencar(1)) / REAL(i2-i1+1))
          grid%z0 = DBLE(SUM(cell(i1:i2)%cencar(2)) / REAL(i2-i1+1))
          IF (debug) WRITE(0,*) 'REMINDER: CHECK GRID CENTROID '//
     .                          'CALCULATION!'
c...        
          CALL BuildConnectionMap(1,nobj)
c...      Mark the radial boundary of the grid - it is assumed that these
c         are relatively simple grids, with well defined nearest 
c         neighbours, i.e. not fully generalised DIVIMP grids that have
c         been tailored:
          DO iobj = 1, nobj
            IF (obj(iobj)%omap(2).EQ.0) obj(iobj)%omap(2) = -1
            IF (obj(iobj)%omap(4).EQ.0) obj(iobj)%omap(4) = -1
          ENDDO

          DEALLOCATE(knot)
          DEALLOCATE(imap)
c          STOP 'whoa'

          CALL osm_DeriveGridQuantities

        CASE (GRD_LOAD_OLD)
          CALL LoadObjects('osm_geometry.raw',status)  ! Change to LoadGeometryObjects...
          IF (status.EQ.-1) 
     .      CALL ER('SetupGrid','Geometry data not found',*99)

          CALL LoadGrid('osm.raw')

          CALL LoadLegacyData('osm_legacy.raw')
        CASE DEFAULT
          CALL ER('ProcessGrid','Unrecognized grid source option',*99)
      ENDSELECT

      CALL DumpData_OSM('output.grid_load','Done loading the grid')

c...
      CALL LoadSupplimentalGridData
      CALL DumpData_OSM('output.grid_sup','Done loading the sup data')

      CALL GenerateTubeGroups
      CALL DumpData_OSM('output.grid_tubes','Done analysing tubes')

      CALL GenerateTargetGroups
      CALL DumpData_OSM('output.grid_targets','Done analysing targets')

      CALL SaveFluidGridGeometry
 
      RETURN
 99   WRITE(0,*) 'IK,IR,I =',ik,ir,i
      STOP
      END
c
c
c
c ======================================================================
c
c
c ======================================================================
c
c subroutine: FindGridCell
c
c     Condition=1 finds two different knots which have the same 
c                 R,Z values for knot number 1. This condition
c                 should only occur at an Xpoint for two cells in the 
c                 inner/outer SOL
c
c     NOTE: The cells sharing vertex 1 or 4 will be in the inner and outer SOL
c           The cells sharing vertex 2 or 3 will be in the core and PFZ
c
c           There seems to be a bug in the code when finding the innermost
c           core ring - using vertex 1 finds a cell in the main SOL adjacent
c           to the core. Stepping inward will usually work if the cell 
c           next to the core is chosen - however, if the separatrix is unusual
c           it could be that the Z coordinate of the cell adjacent to the PFZ
c           will be closer to the center of the plasma and thus the wrong cell
c           will be chosen to find the center of the grid. In order to fix this, 
c           the test vertex for the XPoint should be set to 2 or 3 so that 
c           cells inside and outside the core are chosen - more clearly above 
c           and below the Xpoints. I'll have to check and see if this has an
c           impact on the calculations for double null grids. 
c
c
c     Condition=2 finds knot 1=2 and knot 4=3 for test cell vs. other cells
c
c     Condition=3 finds knot 3=2 and knot 4=1 for test cell vs. other cells
c
c     Condition=4 finds knot 2=1 and knot 3=4 for test cell vs. other cells
c
c     Condition=5 finds knot 1=4 and knot 2=3 for test cell vs. other cells
c
c     Condition=6 finds two different knots which have the same 
c                 R,Z values for knot number 2. This condition
c                 should only occur at an Xpoint for two cells in the core/pfz
c
c
c
      RECURSIVE SUBROUTINE FindGridCell(NUMZONE,izone,
     .                                   condition,index1,index2)
c      RECURSIVE SUBROUTINE FindGridCell(nknot,knot,NUMZONE,izone,
c     .                                   condition,index1,index2)
      USE mod_grid
      IMPLICIT none

      INTEGER index1,index2,NUMZONE,izone(NUMZONE+1,NUMZONE),
     .        condition
c      TYPE(type_grid_cell) :: knot(0:nknot)      

      REAL*8, PARAMETER :: DTOL=1.0D-06

      INTEGER i1,i2,r1,z1
      LOGICAL output

      output = .false.

      i1 = index1

      index2 = -1

      if (condition.ne.1.and.output) then 

         write(6,'(a,100i8)') 'FINDKNOT1:',
     >      nknot,NUMZONE,condition,index1,index2
         write(6,'(a,100i8)') 'FINDKNOT2:',
     >      izone
c
         write(0,'(a,100i8)') 'FINDKNOT1:',
     >      nknot,NUMZONE,condition,index1,index2
c      write(0,'(a,100i8)') 'FINDKNOT:',
c     >      izone

      endif

      DO z1 = knot(i1)%zzone-1, knot(i1)%zzone+1
        IF (z1.LT.1.OR.z1.GT.NUMZONE) CYCLE
        DO r1 = knot(i1)%rzone-1, knot(i1)%rzone+1
          IF (r1.LT.1.OR.r1.GT.NUMZONE) CYCLE  

c          DO i2 = 1, nknot
          DO i2 = izone(r1,z1), izone(r1+1,z1)-1
             
            if (output) then 
               write(6,'(a,15i8)') 'IZONE:',i1,z1,r1,i2,izone(r1,z1),
     >                izone(r1+1,z1)-1
            endif

            IF (i1.EQ.i2) CYCLE
c...
            IF     (condition.EQ.XPT_SEARCH.AND.
     .              DABS(knot(i1)%rv(1)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(1)-knot(i2)%zv(1)).LT.DTOL) THEN

              if (output) then
                 WRITE(0,*) 'XPOINT SOL:',i1,i2
                 WRITE(6,*) 'XPOINT SOL:',i1,i2
              endif
                 
              index2 = i2
              RETURN

            ELSEIF (condition.EQ.LIMITER_SEARCH.AND.
     .              ABS(knot(i1)%ik-knot(i2)%ik).GT.1.AND.
     .              knot(i1)%ir.EQ.knot(i2)%ir.AND.
     .              DABS(knot(i1)%rv(1)-knot(i2)%rv(4)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(1)-knot(i2)%zv(4)).LT.DTOL.AND.
     .              (DABS(knot(i1)%rv(2)-knot(i2)%rv(3)).GT.DTOL.OR.
     .               DABS(knot(i1)%zv(2)-knot(i2)%zv(3)).GT.DTOL)) THEN

c              if (output) then
                 WRITE(0,'(A,6I6)') '------LIMITER POINT:',i1,i2,
     .                   knot(i1)%ik,knot(i1)%ir,
     .                   knot(i2)%ik,knot(i2)%ir
                 WRITE(0,*) knot(i1)%zv(1)
                 WRITE(0,*) knot(i1)%rv(2),knot(i2)%rv(3)
                 WRITE(0,*) knot(i1)%zv(2),knot(i2)%zv(3)
                 WRITE(6,*) 'LIMITER POINT:',i1,i2
c              endif
                 
              index2 = i2
              RETURN

            ELSEIF     (condition.EQ.6.AND.
     .              DABS(knot(i1)%rv(2)-knot(i2)%rv(2)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(2)-knot(i2)%zv(2)).LT.DTOL) THEN
c... *** WHAT'S THE POINT OF THIS? ***
              STOP 'MYSTERIOUS condition.EQ.6'
              if (output) then
                 WRITE(0,*) 'XPOINT CORE/PFZ:',i1,i2
                 WRITE(6,*) 'XPOINT CORE/PFZ:',i1,i2
              endif
                 
              index2 = i2
              RETURN

            ELSEIF (condition.EQ.R_INWARD.AND.
     .              DABS(knot(i1)%rv(1)-knot(i2)%rv(2)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(1)-knot(i2)%zv(2)).LT.DTOL.AND.
     .              DABS(knot(i1)%rv(4)-knot(i2)%rv(3)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(4)-knot(i2)%zv(3)).LT.DTOL) THEN

              if (output) then 
                WRITE(0,*) 'SIDE INWARD 41:',i1,i2
                WRITE(6,*) 'SIDE INWARD 41:',i1,i2
              endif

              index2 = i2
              RETURN

            ELSEIF (condition.EQ.P_FORWARD.AND.
     .              DABS(knot(i1)%rv(3)-knot(i2)%rv(2)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(3)-knot(i2)%zv(2)).LT.DTOL.AND.
     .              DABS(knot(i1)%rv(4)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(4)-knot(i2)%zv(1)).LT.DTOL) THEN
              IF (output) THEN
                WRITE(0,*) 'SIDE UP 34:',i1,i2
                WRITE(6,*) 'SIDE UP 34:',i1,i2
              ENDIF
              index2 = i2
              RETURN
            ELSEIF (condition.EQ.R_OUTWARD.AND.
     .              DABS(knot(i1)%rv(2)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(2)-knot(i2)%zv(1)).LT.DTOL.AND.
     .              DABS(knot(i1)%rv(3)-knot(i2)%rv(4)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(3)-knot(i2)%zv(4)).LT.DTOL) THEN


               if (output) then 
                 WRITE(0,*) 'SIDE OUTWARD 23:',i1,i2
                 WRITE(6,*) 'SIDE OUTWARD 23:',i1,i2
               endif

              index2 = i2
              RETURN

            ELSEIF (condition.EQ.P_BACKWARD.AND.
     .              DABS(knot(i1)%rv(1)-knot(i2)%rv(4)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(1)-knot(i2)%zv(4)).LT.DTOL.AND.
     .              DABS(knot(i1)%rv(2)-knot(i2)%rv(3)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(2)-knot(i2)%zv(3)).LT.DTOL) THEN
c...          Matching sides 12 and 34:

               if (output) then
                 WRITE(0,*) 'SIDE DOWN 12:',i1,i2
                 WRITE(6,*) 'SIDE DOWN 12:',i1,i2
               endif

              index2 = i2
              RETURN



            ENDIF

          ENDDO

        ENDDO
      ENDDO

c
c     If the code reaches here it has failed to find a match in the zoning system - try scanning all cells
c



      RETURN
 99   STOP
      END

c
c ========================================================================
c
      SUBROUTINE MoveGridCell(knot1,knot2)
      USE mod_grid
      IMPLICIT none

c      TYPE type_cell
c        INTEGER :: index,ik,ir,nv,rzone,zzone,xpt,map
c        REAL    :: rcen,zcen,bratio,rv(4),zv(4)
c      ENDTYPE type_cell

      TYPE(type_grid_cell) :: knot1,knot2      

      INTEGER i1

      knot2%index  = knot1%index
      knot2%ik     = knot1%ik
      knot2%ir     = knot1%ir
      knot2%rzone  = knot1%rzone
      knot2%zzone  = knot1%zzone
      knot2%xpt    = knot1%xpt
      knot2%rcen   = knot1%rcen
      knot2%zcen   = knot1%zcen
      knot2%bratio = knot1%bratio
      DO i1 = 1, 4
        knot2%rv(i1) = knot1%rv(i1)
        knot2%zv(i1) = knot1%zv(i1)
      ENDDO

      RETURN
 99   STOP
      END
c     
c     ======================================================================
c     
c     subroutine: LoadGeneralisedGrid
c
c...  (This is not the most efficient way of doing things -- certaintly it would
c     be better just to avoid storing the cells as the grid file is read in --    <---WHAT?
c     but it is general, which is the goal here, and makes minimal assumptions 
c     about the structure of the grid file):
c     

      SUBROUTINE LoadGeneralisedGrid
      USE mod_interface
      USE mod_grid
      USE mod_sol28_global
      IMPLICIT none

c...  Input:
      INTEGER grdfp

c...  Output:
      INTEGER :: irsep,irsep2,irwall,irtrap,nrs,ikti,ikto,nks(1000)
      LOGICAL connected
      
      INTEGER, PARAMETER :: NUMZONE = 5
      REAL*8,  PARAMETER :: HI   = 1.0E+20
      REAL*8,  PARAMETER :: DTOL = 1.0D-07

      INTEGER   i1,i2,z1,r1,kind,nxpt,ixpt(0:2,2),cxpt(0:2,2),i3,
     .          i4,izone(NUMZONE+1,NUMZONE),newi1,icore(0:2,2),id,
     .          tmpnks,istart,fp,maxik,maxir,outfp,count,
     .          ik,ir,irstart,ir1,ir2,idum1,ir_del,nlim,iknot
      LOGICAL   cont,deleteknot,debug,swap,cell_deletion
      REAL*8    vrmin,vzmin,vrmax,vzmax,rspan,zspan,area
      REAL*8    rvdp(4),zvdp(4),areadp
      CHARACTER buffer*1000

      debug = .FALSE.
      outfp  = 0

c...  Read the knot data:

      grd_format   = opt%f_grid_format      ! 1
      grd_filename = TRIM(opt%f_grid_file)  ! 'iterm.carre.105'  ! 'sonnet_13018_250.sm'   
      grdfp = 99
      OPEN(UNIT=grdfp,FILE=TRIM(grd_filename),ACCESS='SEQUENTIAL',
     .     ERR=98)     

      IF (debug) WRITE(0,*) 'grid file name =',TRIM(grd_filename)

      SELECTCASE (grd_format)
c       --------------------------------------------------------------
        CASE (-2:-1)
c...      Find the start of the cell/knot information in the grid file:
          WRITE(buffer,'(1000X)')
          DO WHILE (buffer(4:8).NE.'=====')
             READ(grdfp,'(A10)',END=98) buffer
          ENDDO
c...      Scan the file to see how many cells are in the grid:
          nknot = 0
          maxik = 0
          maxir = 0
          DO WHILE(nknot.EQ.0.OR.buffer(1:7 ).EQ.'Element')
            READ(grdfp,70,END=97) iknot,ik,ir
            READ(grdfp,* ,END=97) 
            READ(grdfp,* ,END=97) 
            READ(grdfp,* ,END=97)
            READ(grdfp,'(A10)',END= 9) buffer
            BACKSPACE(grdfp)
            nknot = nknot + 1
            maxik = MAX(ik+1,maxik)
            maxir = MAX(ir+1,maxir)
          ENDDO
  9       CONTINUE  ! EOF
          nknot = nknot + 1
c...      Setup arrays:
          WRITE(0,*) 'MAXIK,IR',maxik,maxir
          IF (iknot.NE.nknot) THEN
            WRITE(0,*) '***& WARNING *** : NKNOT,IKNOT=',nknot,iknot
          ENDIF
          ALLOCATE(knot(0:nknot))
          ALLOCATE(imap(maxik,0:3*maxir))
          nknot = 0
c...      Load grid:
          REWIND(grdfp)
          WRITE(buffer,'(1000X)')
          DO WHILE (buffer(4:8).NE.'=====')
             READ(grdfp,'(A10)',END=98) buffer
          ENDDO
          DO WHILE(nknot.EQ.0.OR.buffer(1:7 ).EQ.'Element')
c            READ(grdfp,'(A50)',END=19) buffer
c            WRITE(0,*) 'BUFFER:',buffer(1:50)
c            BACKSPACE(grdfp)
            nknot = nknot + 1
            READ(grdfp,70,END=97) knot(nknot)%index,
     .                            knot(nknot)%ik   ,knot(nknot)%ir, 
     .                            knot(nknot)%rv(2),knot(nknot)%zv(2),
     .                            knot(nknot)%rv(3),knot(nknot)%zv(3)
            READ(grdfp,71,END=97) idum1,
     .                            knot(nknot)%rcen ,knot(nknot)%zcen
            READ(grdfp,72,END=97) knot(nknot)%rv(1),knot(nknot)%zv(1),
     .                            knot(nknot)%rv(4),knot(nknot)%zv(4)
            knot(nknot)%bratio = 1.0
          
 70         FORMAT(10X,I5,4X,I6,2x,I6,4x,E16.10,4X,E17.10,7X,E16.10,
     .             4X,E17.10)
 71         FORMAT(13X,E17.10,36X,E16.10,4X,E17.10)
c 71         FORMAT(18X,I1,35X,E16.10,4X,E17.10)
 72         FORMAT(37X,E16.10,4X,E17.10,7X,E16.10,4X,E17.10)
          
c            WRITE(0,*) knot(nknot)%index
c            WRITE(0,*) knot(nknot)%ik,knot(nknot)%ir
c            WRITE(0,*) knot(nknot)%rv(2),knot(nknot)%zv(2)
c            WRITE(0,*) knot(nknot)%rv(3),knot(nknot)%zv(3)
c            WRITE(0,*) knot(nknot)%bratio
c            WRITE(0,*) knot(nknot)%rcen ,knot(nknot)%zcen
c            WRITE(0,*) knot(nknot)%rv(1),knot(nknot)%zv(1)
c            WRITE(0,*) knot(nknot)%rv(4),knot(nknot)%zv(4)
c            STOP 'FSAFDSD'
          
            knot(nknot)%nv = 4
c...        Dividing line in grid file:       
            READ(grdfp,*)
            READ(grdfp,'(A10)',END=19) buffer
            BACKSPACE(grdfp)
          ENDDO
 19       CONTINUE  ! EOF

c       --------------------------------------------------------------
        CASE (GRD_FORMAT_SONNET)
c...      Find the start of the cell/knot information in the grid file:
          WRITE(buffer,'(1000X)')
          DO WHILE (buffer(4:8).NE.'=====')
             READ(grdfp,'(A10)',END=98) buffer
          ENDDO
c...      Scan the file to see how many cells are in the grid:
          nknot = 0
          maxik = 0
          maxir = 0
          DO WHILE(nknot.EQ.0.OR.buffer(4:10).EQ.'Element')
            READ(grdfp,80,END=97) nknot,ik,ir
            READ(grdfp,* ,END=97) 
            READ(grdfp,* ,END=97) 
            READ(grdfp,* ,END=97)
            READ(grdfp,'(A10)',END=10) buffer
            BACKSPACE(grdfp)
            maxik = MAX(ik+1,maxik)
            maxir = MAX(ir+1,maxir)
          ENDDO
 10       CONTINUE  ! EOF
          nknot = nknot + 1
c...      Setup arrays:
          ALLOCATE(knot(0:nknot))
          ALLOCATE(imap(maxik,0:3*maxir))
          nknot = 0
c...      Load grid:
          REWIND(grdfp)
          WRITE(buffer,'(1000X)')
          DO WHILE (buffer(4:8).NE.'=====')
             READ(grdfp,'(A10)',END=98) buffer
          ENDDO
          DO WHILE(nknot.EQ.0.OR.buffer(4:10).EQ.'Element')
c            READ(grdfp,'(A50)',END=10) buffer
c            WRITE(0,*) 'BUFFER:',buffer(1:50)
c            BACKSPACE(grdfp)
            nknot = nknot + 1
            READ(grdfp,80,END=97) knot(nknot)%index,
     .                            knot(nknot)%ik   ,knot(nknot)%ir, 
     .                            knot(nknot)%rv(2),knot(nknot)%zv(2),
     .                            knot(nknot)%rv(3),knot(nknot)%zv(3)
            READ(grdfp,81,END=97) knot(nknot)%bratio,
     .                            knot(nknot)%rcen ,knot(nknot)%zcen
            READ(grdfp,82,END=97) knot(nknot)%rv(1),knot(nknot)%zv(1),
     .                            knot(nknot)%rv(4),knot(nknot)%zv(4)
            knot(nknot)%nv = 4
c...        Dividing line in grid file:       
            READ(grdfp,*)
            READ(grdfp,'(A10)',END=20) buffer
            BACKSPACE(grdfp)
          ENDDO
 80       FORMAT(10X,I5,4X,I3,1x,I3,4x,E17.10,1X,E17.10,8X,E17.10,
     .           1X,E17.10)
 81       FORMAT(18X,E17.10,14x,E17.10,1x,E17.10)
 82       FORMAT(30X,E17.10,1X,E17.10,8x,E17.10,1X,E17.10)
 20       CONTINUE  ! EOF
c       --------------------------------------------------------------
        CASE DEFAULT
          CALL ER('LoadGeneralisedGrid','Unknown grid type',*99)
      ENDSELECT

      IF (debug) THEN
        WRITE(outfp,*) 'GRID LOADED'
        WRITE(outfp,*) '  NKNOT=',nknot
        WRITE(outfp,*) '  MAXIK=',maxik
        WRITE(outfp,*) '  MAXIR=',maxir
      ENDIF

      CALL inOpenInterface('osm.idl.fluid_grid_debug')
      DO i1 = 1, nknot
        CALL inPutData(knot(i1)%ik,'IK','none')
        CALL inPutData(knot(i1)%ir,'IR','none')
        CALL inPutData(4          ,'NV','none')
        DO i2 = 1, 4
          WRITE(buffer,'(A,I1,500X)') 'RV_',i2
          CALL inPutData(knot(i1)%rv(i2),TRIM(buffer),'m')
          WRITE(buffer,'(A,I1,500X)') 'ZV_',i2
          CALL inPutData(knot(i1)%zv(i2),TRIM(buffer),'m')
        ENDDO
      ENDDO
      CALL inCloseInterface
 
c...  Delete zero volume cells:
c     
      IF (debug) WRITE(outfp,*) 'REMOVING ZERO VOLUME CELLS'
      cell_deletion = .FALSE.
      i1 = 1
      DO WHILE(i1.LE.nknot)
        IF     (opt%f_grid_strip.EQ.1.AND.
     .          (knot(i1)%ik.EQ.0.OR.knot(i1)%ik.EQ.maxik-1.OR. ! virtual cells on end of rings
     .           knot(i1)%ir.EQ.0.OR.knot(i1)%ir.EQ.maxir-1)) THEN ! virtual rings
c...      Forced removal of boundary cells (nuclear option):
          deleteknot = .TRUE.
        ELSEIF ((DABS(knot(i1)%rv(3)-knot(i1)%rv(2)).LT.DTOL.AND.
     .           DABS(knot(i1)%zv(3)-knot(i1)%zv(2)).LT.DTOL.AND.
     .           DABS(knot(i1)%rv(4)-knot(i1)%rv(1)).LT.DTOL.AND.
     .           DABS(knot(i1)%zv(4)-knot(i1)%zv(1)).LT.DTOL).OR.
     .          (DABS(knot(i1)%rv(2)-knot(i1)%rv(1)).LT.DTOL.AND.
     .           DABS(knot(i1)%zv(2)-knot(i1)%zv(1)).LT.DTOL.AND.
     .           DABS(knot(i1)%rv(4)-knot(i1)%rv(3)).LT.DTOL.AND.
     .           DABS(knot(i1)%zv(4)-knot(i1)%zv(3)).LT.DTOL)) THEN
c...      Also get rid of zero volume cells, which can be present in UEDGE
c         double null grids.  The above condition is the best identifier
c         for these (for grids generated with UEDGE anyway):
          deleteknot = .TRUE.
        ELSE
c...      Cell to be kept, advance index:
          deleteknot = .FALSE.
          i1 = i1 + 1
        ENDIF
c...    The request has been made to delete the current cell:
        IF (deleteknot) THEN
          cell_deletion = .TRUE.
          IF (i1.LT.nknot) THEN
            DO i2 = i1, nknot-1
              CALL MoveGridCell(knot(i2+1),knot(i2))
            ENDDO
          ENDIF
          nknot = nknot - 1
        ENDIF
      ENDDO
      IF (debug) THEN
        IF (cell_deletion) THEN
          WRITE(outfp,*) '  NKNOT=',nknot
        ELSE
          WRITE(outfp,*) '  No cells deleted'
        ENDIF
      ENDIF
c...  Cell integrity checks (need to come up with something better -- like
c     something that looks for misshapen cells, i.e. XBC's):
      DO i1 = 1, nknot
        IF ((DABS(knot(i1)%rv(3)-knot(i1)%rv(2)).LT.DTOL.AND.
     .       DABS(knot(i1)%zv(3)-knot(i1)%zv(2)).LT.DTOL).OR.
     .      (DABS(knot(i1)%rv(4)-knot(i1)%rv(1)).LT.DTOL.AND.
     .       DABS(knot(i1)%zv(4)-knot(i1)%zv(1)).LT.DTOL).OR.
     .      (DABS(knot(i1)%rv(3)-knot(i1)%rv(1)).LT.DTOL.AND.
     .       DABS(knot(i1)%zv(3)-knot(i1)%zv(1)).LT.DTOL).OR.
     .      (DABS(knot(i1)%rv(4)-knot(i1)%rv(2)).LT.DTOL.AND.
     .       DABS(knot(i1)%zv(4)-knot(i1)%zv(2)).LT.DTOL)) THEN
c...      Triangles, maybe:
          CALL ER('LoadGeneralisedGrid','Possible triangular cell '//
     .            'detected',*99)
        ENDIF          
      ENDDO

c...  Global R,Z shifts:
c       Not implemented yet...

c...  Assign cell zone assignment based on its spatial location in the
c     grid.  This (dramatically) increases the speed of the search
c     algorithm when building the connection map (grid structure) below.
      vrmin =  HI
      vrmax = -HI
      vzmin =  HI
      vzmax = -HI
c...  Find the spatial extent of the grid:
      DO i1 = 1, nknot      
         DO i2 = 1, knot(i1)%nv
            IF (knot(i1)%rv(i2).LT.vrmin) vrmin = knot(i1)%rv(i2)
            IF (knot(i1)%rv(i2).GT.vrmax) vrmax = knot(i1)%rv(i2)
            IF (knot(i1)%zv(i2).LT.vzmin) vzmin = knot(i1)%zv(i2)
            IF (knot(i1)%zv(i2).GT.vzmax) vzmax = knot(i1)%zv(i2)
         ENDDO
      ENDDO
      vrmin = vrmin - 0.001D0  ! Expand the domain slightly...
      vrmax = vrmax + 0.001D0
      vzmin = vzmin - 0.001D0
      vzmax = vzmax + 0.001D0
      IF (debug) THEN
        WRITE(outfp,*) 'GRID DOMAIN EXTENT:'
        WRITE(outfp,*) '  ',vrmin
        WRITE(outfp,*) '  ',vrmax
        WRITE(outfp,*) '  ',vzmin
        WRITE(outfp,*) '  ',vzmax
      ENDIF
c...  Assign each cell to a zone:
      rspan = (vrmax - vrmin) / DBLE(NUMZONE)
      zspan = (vzmax - vzmin) / DBLE(NUMZONE)
      DO i1 = 1, nknot
        knot(i1)%rzone = INT( (knot(i1)%rcen - vrmin) / rspan ) + 1
        knot(i1)%zzone = INT( (knot(i1)%zcen - vzmin) / zspan ) + 1
      ENDDO
c...  Sort the cells so that they are grouped according to which zone:
      kind = 1
      DO z1 = 1, NUMZONE
        DO r1 = 1, NUMZONE
          izone(r1,z1) = kind
          DO i1 = kind, nknot
            IF (knot(i1)%rzone.EQ.r1.AND.knot(i1)%zzone.EQ.z1) THEN
              IF (i1.EQ.kind) THEN
c...            Do nothing:
              ELSE
c...            Swap the cells:
                CALL MoveGridCell(knot(kind),knot(0)   )
                CALL MoveGridCell(knot(i1)  ,knot(kind))
                CALL MoveGridCell(knot(0)   ,knot(i1)  )
              ENDIF
              kind = kind + 1
            ENDIF
          ENDDO
        ENDDO
      ENDDO
c...  
      DO i1 = 1, NUMZONE-1
         izone(NUMZONE+1,i1) = izone(1,i1+1)
      ENDDO
      izone(NUMZONE+1,NUMZONE) = nknot + 1
      IF (debug) THEN
        WRITE(outfp,*) 'ZONE POPULATIONS:'
        DO i2 = 1, NUMZONE
          WRITE(0,'(4X,I6,2X,10I6)') i2,izone(1:NUMZONE+1,i2)
        ENDDO
      ENDIF
c...  

c...  Find and order the x-point cells:
c     -------------------------------------------------------------------     
      nxpt = 0
      DO i1 = 1, nknot
        IF (knot(i1)%xpt.NE.0) CYCLE
        CALL FindGridCell(NUMZONE,izone,XPT_SEARCH,i1,i2)
        IF (i2.NE.-1) THEN
          nxpt = nxpt + 1
          ixpt(nxpt,1) = i1
          ixpt(nxpt,2) = i2
          knot(i1)%xpt = i2
          knot(i2)%xpt = i1
        ENDIF
      ENDDO
c...  No x-points found, for some reason:
      IF (debug) THEN
        WRITE(outfp,*) 'X-POINT CELLS:'
        WRITE(outfp,*) '  IXPT1=',ixpt(1:nxpt,1)
        WRITE(outfp,*) '  IXPT2=',ixpt(1:nxpt,2)
      ENDIF
c...  Check if it's a limiter grid:          
      nlim = 0
      IF (nxpt.EQ.0) THEN 
        DO i1 = 1, nknot
          CALL FindGridCell(NUMZONE,izone,LIMITER_SEARCH,i1,i2)
          IF (i2.NE.-1) THEN
            nlim = nlim + 1
            ixpt(nlim,1) = i1
            ixpt(nlim,2) = i2
          ENDIF
        ENDDO        
      ENDIF

      IF (nxpt.EQ.0.AND.nlim.EQ.0) 
     .  CALL ER('LoadGeneralisedgrid','No x- or limiter points '//
     .          'found',*99)
      IF (nxpt.GT.2) 
     .  CALL ER('LoadGeneralizedGrid','More than 2 x-points '//
     .          'identified, not good',*99)
      IF (nlim.GT.1) 
     .  CALL ER('LoadGeneralizedGrid','More than 1 limiter points '//
     .          'identified, not good',*99)

c...  
      cxpt = 0
      icore = 0
      DO i4 = 1, MAX(nxpt,nlim)
        DO i3 = 1, 2
          newi1 = ixpt(i4,i3)
          cont = .TRUE.
          DO WHILE (cont)
            i1 = newi1
            cxpt(i4,i3) = cxpt(i4,i3) + 1
            cont = .FALSE.
            CALL FindGridCell(NUMZONE,izone,R_INWARD,i1,i2)
            IF (i2.NE.-1) THEN
               cont = .TRUE.
               newi1 = i2 
               icore(i4,i3) = i2
            ENDIF
          ENDDO
c...      Now, somewhat labouriously, need to find out if the current focus
c         cell ends on a continuous ring -- if yes, then this is the 
c         principal cell for this x-point pair, i.e. the one adjacent to the core:
          i1 = icore(i4,i3)
          count = 0
          DO WHILE (count.LT.5000)  
            count = count + 1
            CALL FindGridCell(NUMZONE,izone,P_FORWARD,i1,i2)  ! Move forward on a ring
            IF     (i2.EQ.-1) THEN
c             Neighbouring cell not found, so at the end of a ring, hence 
c             somewhere in a private flux region:
              icore(i4,i3) = -icore(i4,i3)  ! Switch sign to mark the result
              EXIT
            ELSEIF (i2.EQ.icore(i4,i3)) THEN
c             Have returned to the starting cell, so the ring must be closed, 
c             i.e. be in the core:
              EXIT
            ENDIF
            i1 = i2
          ENDDO
          IF (count.EQ.5000) 
     .       CALL ER('LoadGeneralisedGrid','Confusing result when '//
     .               'processing x/limiter-points (bad grid?)',*99)
        ENDDO
c...    Swap the cells, if necessary:
        IF (nxpt.GT.0.AND.icore(i4,1).LT.0) THEN
          ixpt (0 ,1) = ixpt (i4,1)
          ixpt (i4,1) = ixpt (i4,2)
          ixpt (i4,2) = ixpt (0 ,1)
          cxpt (0 ,1) = cxpt (i4,1)
          cxpt (i4,1) = cxpt (i4,2)
          cxpt (i4,2) = cxpt (0 ,1)
          icore(0 ,1) = icore(i4,1)
          icore(i4,1) = icore(i4,2)
          icore(i4,2) = icore(0 ,1)
        ENDIF
      ENDDO
      IF (debug) THEN
        WRITE(outfp,*) 'X/LIMITER-POINT SCANNING:' 
        WRITE(outfp,*) '   IXTP1 =',ixpt (1:MAX(nxpt,nlim),1)
        WRITE(outfp,*) '   CXPT1 =',cxpt (1:MAX(nxpt,nlim),1)
        WRITE(outfp,*) '   ICORE1=',icore(1:MAX(nxpt,nlim),1)
        WRITE(outfp,*) 
        WRITE(outfp,*) '   IXTP2 =',ixpt (1:MAX(nxpt,nlim),2)
        WRITE(outfp,*) '   CXPT2 =',cxpt (1:MAX(nxpt,nlim),2)
        WRITE(outfp,*) '   ICORE2=',icore(1:MAX(nxpt,nlim),2)
      ENDIF
c...  Check that the x-points are ordered properly, with the primary x-point
c     at index 1, and whether or not the double-null grid is connected:
      connected = .FALSE.
      IF (nxpt.GT.1) THEN
        swap = .FALSE.
        IF     (nxpt.GT.1.AND.cxpt(1,1).EQ.cxpt(2,1)) THEN
c...      Connected double-null grid indicated -- want lower x-point to 
c         be "primary", as a convention (the "separatrix" in the code runs 
c         up the HFS (LHS) of the core:
          IF (knot(ixpt(1,1))%zcen.GT.knot(ixpt(2,1))%zcen) swap=.TRUE.  ! 
          connected = .TRUE.
        ELSEIF (nxpt.GT.1.AND.cxpt(1,1).GT.cxpt(2,1)) THEN
          swap = .TRUE.
        ENDIF
        IF (swap) THEN
          ixpt (0,1:2) = ixpt (1,1:2)
          ixpt (1,1:2) = ixpt (2,1:2)
          ixpt (2,1:2) = ixpt (0,1:2)
          cxpt (0,1:2) = cxpt (1,1:2)
          cxpt (1,1:2) = cxpt (2,1:2)
          cxpt (2,1:2) = cxpt (0,1:2)
          icore(0,1:2) = icore(1,1:2)
          icore(1,1:2) = icore(2,1:2)
          icore(2,1:2) = icore(0,1:2)
          IF (debug) THEN
            WRITE(outfp,*) 'X/LIMITER-POINT SWAP:' 
            WRITE(outfp,*) '   IXTP1 =',ixpt (1:MAX(nxpt,nlim),1)
            WRITE(outfp,*) '   CXPT1 =',cxpt (1:MAX(nxpt,nlim),1)
            WRITE(outfp,*) '   ICORE1=',icore(1:MAX(nxpt,nlim),1)
            WRITE(outfp,*) 
            WRITE(outfp,*) '   IXTP2 =',ixpt (1:MAX(nxpt,nlim),2)
            WRITE(outfp,*) '   CXPT2 =',cxpt (1:MAX(nxpt,nlim),2)
            WRITE(outfp,*) '   ICORE2=',icore(1:MAX(nxpt,nlim),2)
          ENDIF
        ENDIF
      ENDIF

c...  Assemble the grid:
c     ------------------------------------------------------------------

c...  Location of the primary separatrix is known:
      irsep  = cxpt(1,1)
      irsep2 = irsep

c...  Build the core rings:
c     ------------------------------------------------------------------
      IF (debug) WRITE(outfp,*) 'PROCESSING CORE RINGS'
      nrs = 1
      ik = 1
      ir = 1
      i1 = icore(1,1) 
      DO WHILE(ir.LT.irsep)
        imap(1,ir) = i1
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
          CALL FindGridCell(NUMZONE,izone,P_FORWARD,i1,i2)  ! Move forward
          IF (i2.NE.-1) THEN
            IF (i2.NE.imap(1,ir)) THEN
              i1 = i2 
              ik = ik + 1
              imap(ik,ir) = i1
              cont = .TRUE.
c             IF (debug) WRITE(0,*) 'CORE MAP:',ik,ir,i1
            ELSE
c             Ring has closed on itself, finished here:
            ENDIF
          ELSE
            CALL ER('LoadGeneralisedgrid','Cell sequence broken in '//
     .              'the core',*99)
          ENDIF
        ENDDO
        nks(ir) = ik
c...    Step outward to the next ring, still in the core:        
        CALL FindGridCell(NUMZONE,izone,R_OUTWARD,imap(1,ir),i2)  ! And to the left
        IF (i2.NE.-1) THEN        
          i1 = i2
          ik = 1
          ir = ir + 1
        ELSE
          CALL ER('LoadGeneralisedgrid','Radial step has failed',*99)
        ENDIF
      ENDDO
      nrs = ir
      IF (debug) THEN
        WRITE(outfp,*) ' CORE RINGS:'
        WRITE(outfp,*) '  NRS  =',nrs
        WRITE(outfp,*) '  IRSEP=',irsep
        WRITE(outfp,*) '  NKS  =',nks(1:irsep-1)
      ENDIF

c...  SOL rings:
c     ------------------------------------------------------------------
      IF (debug) WRITE(outfp,*) 'PROCESSING SOL RINGS'
c...  Step out of the core:
      i1 = imap(1,irsep-1)
      CALL FindGridCell(NUMZONE,izone,R_OUTWARD,i1,i2)
      IF (i2.NE.-1) THEN  
c...    We should be back at the primary x-point:
        IF (i2.NE.ixpt(1,1)) 
     .    CALL ER('LoadGeneralisedGrid','Not back at the x/limiter-'//
     .            'point cell for some reason',*99)
c...    Move backward along the ring to find the target cell:
        i1 = i2
        DO WHILE(i2.NE.-1)
          CALL FindGridCell(NUMZONE,izone,P_BACKWARD,i1,i2) ! Moving backward 
          IF (i2.NE.-1) i1 = i2
        ENDDO
      ELSE
        CALL ER('LoadGeneralisedGrid','Radial step out from the '//
     .          'core failed',*99)
      ENDIF
c...  Target located, start mapping the SOL:
      ik = 1
      ir = irsep
      imap(ik,ir) = i1
      cont = .TRUE.
      DO WHILE(cont)
        cont = .FALSE.
c...    Move forward along the ring until the corresponding target is reached:
        CALL FindGridCell(NUMZONE,izone,P_FORWARD,i1,i2)
        IF (i2.NE.-1) THEN
          i1 = i2 
          ik = ik + 1
          imap(ik,ir) = i1
          cont = .TRUE.
c          IF (debug) WRITE(0,*) 'SOL MAP (INNER):',ik,ir,i1
        ENDIF
        IF (.NOT.cont) THEN
c...      Step radially outward if ring is finished, i.e. a cell has been reached
c         for which no "next" cell is found when moving along the ring:
          nks(ir) = ik
          i1 = imap(1,ir)
          CALL FindGridCell(NUMZONE,izone,R_OUTWARD,i1,i2)  ! Step radially          
          IF (i2.NE.-1) THEN
             i1 = i2
             ik = 1
             ir = ir + 1
             imap(ik,ir) = i1
             cont = .TRUE.
          ELSE
c...        Finished here, no outward radial setup found:
            EXIT 
c            CALL ER('LoadGeneralisedgrid','Radial step '//
c     .              'failed in the SOL',*99)
          ENDIF
        ENDIF
      ENDDO
      irwall = ir
      irtrap = ir + 1
      nrs = ir
      IF (debug) THEN
        WRITE(outfp,*) ' SOL RINGS:'
        WRITE(outfp,*) '  IRSEP =',irsep
        WRITE(outfp,*) '  IRWALL=',irwall
        WRITE(outfp,*) '  IRTRAP=',irtrap
        WRITE(outfp,*) '  NRS   =',nrs
      ENDIF
c...  SOL rings associated with double-null grids, including the 
c     secondary private flux region:
c     ------------------------------------------------------------------
      IF (nxpt.GT.1) THEN
c...    Find scrape-off layer rings that were not processed yet, which 
c       happens with double-null grids:

        irsep2 = irsep + cxpt(2,1) - cxpt(1,1) - 1

        IF (debug) WRITE(outfp,*) 'PROCESSING SECONDARY PFR'
c...    Process the secondary x-point PFR, which is just considered part of the
c       SOL for generalized grids:
        i1 = ixpt(2,2)
c...    Step radially from the secondary x-point into the PFR:
        CALL FindGridCell(NUMZONE,izone,R_INWARD,i1,i2)
        IF (i2.EQ.-1) CALL ER('LoadGeneralisedgrid','Unable to step '//
     .                        'from the secondary x-pt to the PFR',*99)
        i1 = i2
        ik = 1
        ir = nrs + 1
        imap(ik,ir) = i1
        irstart = ir
c...    Proceed alon the ring to the target:
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
          CALL FindGridCell(NUMZONE,izone,P_BACKWARD,i1,i2)
          IF (i2.NE.-1) THEN 
            i1 = i2
            ik = ik + 1
            imap(ik,ir) = i1
            cont = .TRUE.
          ENDIF
        ENDDO
c...    Start mapping:
        ik = 1
        ir = nrs + 1
        imap(ik,ir) = i1
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
          CALL FindGridCell(NUMZONE,izone,P_FORWARD,i1,i2)
          IF (i2.NE.-1) THEN
            i1 = i2 
            ik = ik + 1
            imap(ik,ir) = i1
            cont = .TRUE.
          ENDIF
c...      Step radially inward if ring is finished:
          IF (.NOT.cont) THEN
            nks(ir) = ik
            i1 = imap(1,ir)
            CALL FindGridCell(NUMZONE,izone,R_INWARD,i1,i2)          
            IF (i2.NE.-1) THEN
              i1 = i2
              ik = 1
              ir = ir + 1
              imap(ik,ir) = i1
              cont = .TRUE.
            ELSE
c...          Assume the outer boundary of the grid:
            ENDIF
          ENDIF
        ENDDO
        irwall = ir
        irtrap = ir + 1
        nrs    = ir
c...    Register assigned knots:
        DO ir = 1, nrs
          knot(imap(1:nks(ir),ir))%map = 1
        ENDDO
c...    Need to reorder the rings:
        DO i1 = 0, (nrs-irstart+1)/2-1
          ir1 = irstart + i1
          ir2 = nrs     - i1
          tmpnks   = nks(ir1)
          nks(ir1) = nks(ir2)
          nks(ir2) = tmpnks
          imap(:,0  ) = imap(:,ir1)
          imap(:,ir1) = imap(:,ir2)
          imap(:,ir2) = imap(:,0  )
        ENDDO

        IF (debug) WRITE(outfp,*) 'DOUBLE NULL GRID, PROCESSING '// 
     .                             'REMAINDER OF SOL'

c...    Start from the last cell on the last core ring:
        i1=imap(nks(irsep-1),irsep-1) 

c...    Step radially outward into the SOL:
        CALL FindGridCell(NUMZONE,izone,R_OUTWARD,i1,i2)
        IF (i2.NE.-1) THEN
          i1 = i2
          IF (knot(i2)%map.EQ.1) THEN  
c...        Keep stepping outward until a cell with no assigned 
c           mapping is found:
            cont = .TRUE.
            DO WHILE(cont)
              cont = .FALSE.
              CALL FindGridCell(NUMZONE,izone,R_OUTWARD,i1,i2)
              IF (i2.NE.-1) THEN 
                i1 = i2
                IF (knot(i1)%map.NE.0) cont = .TRUE.
              ELSE
                CALL ER('LoadGeneralisedGrid','Unable to locate '//
     .                  'unassigned cell in the SOL',*99)
              ENDIF
            ENDDO
          ENDIF
        ELSE
          CALL ER('LoadGeneralisedGrid','Bad radial step from '//
     .            'core to SOL',*99)
        ENDIF
c...    Proceed backward along the ring to the target:
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
          CALL FindGridCell(NUMZONE,izone,P_BACKWARD,i1,i2)
          IF (i2.NE.-1) THEN 
            i1 = i2
            cont = .TRUE.
          ENDIF
        ENDDO
c...    Start mapping:
        ik = 1
        ir = nrs + 1
        IF (connected) irsep2 = ir
        imap(ik,ir) = i1
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
c...      Move forward along the ring:
          CALL FindGridCell(NUMZONE,izone,P_FORWARD,i1,i2)
          IF (i2.NE.-1) THEN
            i1 = i2 
            ik = ik + 1
            imap(ik,ir) = i1
            cont = .TRUE.
          ENDIF
c...      Step radially outward if ring is finished:
          IF (.NOT.cont) THEN
            nks(ir) = ik
            i1 = imap(1,ir)
            CALL FindGridCell(NUMZONE,izone,R_OUTWARD,i1,i2)          
            IF (i2.NE.-1) THEN
              i1 = i2
              ik = 1
              ir = ir + 1
              imap(ik,ir) = i1
              cont = .TRUE.
            ELSE
c...          Assume the outer boundary of the grid:
            ENDIF
          ENDIF
        ENDDO
        irwall = ir
        irtrap = ir + 1
        nrs = ir
c...    Register all knots that have been mapped to the grid, again:
        DO ir = 1, nrs
          knot(imap(1:nks(ir),ir))%map = 1
        ENDDO

      ENDIF  ! Done processing double-null rings
c
c...  Primary private flux zone (PFZ) rings:
c     ------------------------------------------------------------------
      IF (nxpt.GE.1) THEN
        IF (debug) WRITE(outfp,*) 'PROCESSING PRIMARY PFZ'

c...    Process the primary x-point PFR:
        i1 = ixpt(1,2)
c...    Move radially from the SOL into the PFR:
        CALL FindGridCell(NUMZONE,izone,R_INWARD,i1,i2)
        IF (i2.EQ.-1) CALL ER('LoadGeneralisedgrid','Failed to step '//
     .                        'into the PFR from the SOL',*99)
        i1 = i2
c...    Move backward on the ring until the target is found:
        cont  = .TRUE.
        DO WHILE(cont)
          cont  = .FALSE.
          CALL FindGridCell(NUMZONE,izone,P_BACKWARD,i1,i2)
          IF (i2.NE.-1) THEN 
            i1 = i2
            cont = .TRUE.
          ENDIF
        ENDDO
c...    Target located, start mapping the primary PFR:
        ik = 1
        ir = nrs + 1
        imap(ik,ir) = i1
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
c...      Move forward along the ring:
          CALL FindGridCell(NUMZONE,izone,P_FORWARD,i1,i2)
          IF (i2.NE.-1) THEN
            i1 = i2 
            ik = ik + 1
            imap(ik,ir) = i1
            cont = .TRUE.
          ENDIF
c...      Step radially outward when the ring is finished:
          IF (.NOT.cont) THEN
            nks(ir) = ik
            i1 = imap(1,ir)
            CALL FindGridCell(NUMZONE,izone,R_INWARD,i1,i2)          
            IF (i2.NE.-1) THEN
               i1 = i2
               ik = 1
               ir = ir + 1
               imap(ik,ir) = i1
               cont = .TRUE.
            ELSE
c...          Assume the outer boundary of the grid:
            ENDIF
          ENDIF
        ENDDO
        nrs = ir

c...    Need to reorder the rings in the primary PFZ:
        DO i1 = 0, (nrs-irtrap+1)/2-1
          ir1 = irtrap + i1
          ir2 = nrs    - i1
          tmpnks   = nks(ir1)
          nks(ir1) = nks(ir2)
          nks(ir2) = tmpnks
          imap(:,0  ) = imap(:,ir1)
          imap(:,ir1) = imap(:,ir2)
          imap(:,ir2) = imap(:,0  )
        ENDDO
      ENDIF

c...  Find IKTO,IKTI:
      IF (nlim.EQ.1) THEN
        ikto = 0
        ikti = nks(irsep) + 1
      ELSE
        ikto = -1
        ikti = -1
        DO ik = 1, nks(irsep)
          IF (connected) THEN
            IF (imap(ik,irsep).EQ.ixpt(1,1) ) ikto = ik 
            IF (imap(ik,irsep).EQ.ixpt(2,2) ) ikti = ik - 1
          ELSE
            IF (imap(ik,irsep).EQ.ixpt(1,1) ) ikto = ik - 1
            IF (imap(ik,irsep).EQ.ixpt(1,2) ) ikti = ik 
          ENDIF
        ENDDO
        IF (ikto.EQ.-1.OR.ikti.EQ.-1)
     .       CALL ER('LoadGeneralisedgrid','IKTI or IKTO not found',*99)
      ENDIF

c.... Delete some rings if necessary (easy to do here since not much 
c     geometry processing has been done yet):
      DO i1 = 1, opt%grd_ntdel
        ir_del = opt%grd_tdel(i1)
        DO ir = ir_del, nrs-1
          nks(ir) = nks(ir+1)
          imap(1:nks(ir),ir) = imap(1:nks(ir),ir+1)
        ENDDO
        IF (ir_del.LT.irsep ) irsep  = irsep  - 1
        IF (ir_del.LT.irsep2) irsep2 = irsep2 - 1
        IF (ir_del.LT.irwall) irwall = irwall - 1
        IF (ir_del.LT.irtrap) irtrap = irtrap - 1
        nrs = nrs - 1
      ENDDO




c...  Assign data to the global arrays:
      grid_load%irsep      = irsep      
      grid_load%irsep2     = irsep2     
      grid_load%irwall     = irwall     
      grid_load%irtrap     = irtrap     
      grid_load%nrs        = nrs        
      grid_load%ikti       = ikti       
      grid_load%ikto       = ikto       
      grid_load%nks(1:nrs) = nks(1:nrs)





      CALL inOpenInterface('osm.idl.fluid_grid_debug')
      DO ir = 1, nrs
        CALL inPutData(nks(ir),'NKS','NA')    
        DO ik = 1, nks(ir)
          CALL inPutData(ik,'IK','NA')    
          CALL inPutData(ir,'IR','NA')    
          CALL inPutData(4 ,'NV','NA')
          i1 = imap(ik,ir)
          DO i2 = 1, 4
            WRITE(buffer,'(A,I1,500X)') 'RV_',i2
            CALL inPutData(knot(i1)%rv(i2),TRIM(buffer),'m')
            WRITE(buffer,'(A,I1,500X)') 'ZV_',i2
            CALL inPutData(knot(i1)%zv(i2),TRIM(buffer),'m')
          ENDDO
        ENDDO
      ENDDO
      CALL inCloseInterface

      RETURN
 96   CALL ER('LoadGeneralisedgrid','Grid file not found',*99)
 97   CALL ER('LoadGeneralisedgrid','Unexpected end-of-file',*99)
 98   CALL ER('LoadGeneralisedgrid','Problem accessing grid file',*99)
 99   WRITE(0,*) 'GRID: ',TRIM(grd_filename)
      WRITE(0,*) 'IXPT: ',ixpt(1,1),ixpt(2,1)
      WRITE(0,*) 'IR  : ',ir
      STOP  
      END
c


c
c ======================================================================
c
      SUBROUTINE DynamicMap(map_iobj,map_iside,map_icell,map_itube)
      USE mod_geometry
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      INTEGER, INTENT(IN)  :: map_iobj ,map_iside
      INTEGER, INTENT(OUT) :: map_icell,map_itube

      INTEGER GetObject,GetTube

      INTEGER  i1,i2,ic1,ic2,itube,icell,iobj,isrf,ivtx(2)
      REAL*8   a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd,p(3),maxtab

      map_icell = -1
      map_itube = -1

      isrf      = ABS(obj(map_iobj)%iside(map_iside))
      ivtx(1:2) = srf(isrf)%ivtx(1:2)
      CALL CalcCentroid(map_iobj,2,p)
      a1 = p(1)
      a2 = p(2)
      b1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
      b2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2))) 

      maxtab = 1.0D+20
      DO itube = 1, ntube
        IF (GetTube(map_iobj,IND_OBJECT).EQ.itube) CYCLE
        ic1 = tube(itube)%cell_index(LO)
        ic2 = tube(itube)%cell_index(HI)
        DO icell = ic1, ic2
          iobj = GetObject(icell,IND_CELL)
          isrf = ABS(obj(iobj)%iside(1))
          ivtx(1:2) = srf(isrf)%ivtx(1:2)
          c1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
          c2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2)))
          isrf = ABS(obj(iobj)%iside(3))
          ivtx(1:2) = srf(isrf)%ivtx(1:2)
          d1 = 0.5D0 * (vtx(1,ivtx(1)) + vtx(1,ivtx(2)))
          d2 = 0.5D0 * (vtx(2,ivtx(1)) + vtx(2,ivtx(2)))
          CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
          IF (tab.GE.0.0D0.AND.tab.LT.maxtab.AND.
     .        tcd.GE.0.0D0.AND.tcd.LT.1.0D0) THEN
            maxtab = tab
            map_icell = icell
            map_itube = itube
          ENDIF
        ENDDO
      ENDDO

      WRITE(88,*) 'DYNAMIC:'
      WRITE(88,*) '  : ',map_iobj ,map_iside
      WRITE(88,*) '  : ',a1,a2
      WRITE(88,*) '  : ',b1,b2
      WRITE(88,*) '  : ',map_icell,map_itube

      RETURN
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE ClipWallToGrid(nwall,rwall,zwall,MAXNSEG,delete)
      USE mod_geometry
      USE mod_sol28_global
      IMPLICIT none

      INTEGER, INTENT(IN)    :: MAXNSEG
      LOGICAL, INTENT(IN)    :: delete
      INTEGER, INTENT(INOUT) :: nwall
      REAL*8 , INTENT(INOUT) :: rwall(MAXNSEG,2),zwall(MAXNSEG,2)

      INTEGER GetTube       
 
      INTEGER, PARAMETER :: MAXNLIST = 1000
      REAL*8 , PARAMETER :: DTOL = 1.0D-07

      INTEGER fp,iobj,itube,nlist,ilist(MAXNLIST,2),clist(MAXNLIST,2),
     .        tube_set,i1,i2,i3,swall(nwall),iwall,mlist(MAXNLIST)
      LOGICAL debug
      REAL*8  x1,x2,y1,y2,x3,x4,y3,y4,s12,s34,s12max,length,
     .        xlist(MAXNLIST,2),ylist(MAXNLIST,2),store_x2,store_y2

      fp = 88
      debug = .FALSE.

      CALL DumpData_OSM('output.clipping','Trying to clip grid')

c...  Assume wall is clockwise-specified.  Cut the grid by extending 
c     the poloidal surfaces of end cells that lie on the radial 
c     fluid grid boundary surfaces.  Delete everything that's between
c     cuts, as appropriate:

      tube_set = -1
      nlist = 0

      DO iobj = 1, nobj
        IF (grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID.OR.  ! Only for fluid grid objects
     .      grp(obj(iobj)%group)%type  .NE.GRP_QUADRANGLE   .OR.
     .      GetTube(iobj,IND_OBJECT).LT.grid%isep) CYCLE

        IF (obj(iobj)%omap(2).EQ.-1.OR.obj(iobj)%omap(4).EQ.-1) THEN
          itube = GetTube(iobj,IND_OBJECT)
          IF     (tube_set.NE.itube) THEN
            tube_set = itube         ! There's an assumption here that 
            nlist = nlist + 1        ! the cells on a ring are sequential
            ilist(nlist,1:2) = iobj  ! which should be fine... - SL, 01/12/09
          ELSEIF (tube_set.EQ.itube) THEN
            ilist(nlist,2  ) = iobj
          ENDIF
        ENDIF
      ENDDO

      IF (debug) THEN
        DO i1 = 1, nlist
          WRITE(fp,*) 'CLIP LIST:',ilist(i1,1:2),nwall
        ENDDO
      ENDIF

c...  Collect the cuts:
      clist = 0
      mlist = 0
      xlist = 0.0D0
      ylist = 0.0D0
      i1 = 0
      DO WHILE(i1.LT.nlist)
        i1 = i1 + 1
        DO i2 = 1, 2
          iobj = ilist(i1,i2)
          IF (obj(iobj)%omap(2).EQ.-1) THEN
            mlist(i1) = 1
            IF (i2.EQ.1) THEN
              CALL GetVertex(iobj,3,x1,y1)
              CALL GetVertex(iobj,2,x2,y2)
            ELSE
              CALL GetVertex(iobj,2,x1,y1)
              CALL GetVertex(iobj,3,x2,y2)
            ENDIF
          ELSE
            mlist(i1) = 2
            IF (i2.EQ.1) THEN
              CALL GetVertex(iobj,4,x1,y1)
              CALL GetVertex(iobj,1,x2,y2)
            ELSE
              CALL GetVertex(iobj,1,x1,y1)
              CALL GetVertex(iobj,4,x2,y2)
            ENDIF
          ENDIF
c         Increase the length of the line segment in case there's
c         a nominal mismatch between the wall and target 
c         specifications or if the line segment is very short:
          store_x2 = x2
          store_y2 = y2
          length = DSQRT((x1-x2)**2 + (y1-y2)**2)
          x2 = x1 + MAX(2.0D0,0.1D0 / length) * (x2 - x1)
          y2 = y1 + MAX(2.0D0,0.1D0 / length) * (y2 - y1)
 1        x1 = store_x2 + MAX(2.0D0,0.1D0 / length) * (x1 - store_x2)
          y1 = store_y2 + MAX(2.0D0,0.1D0 / length) * (y1 - store_y2)

          IF (debug) THEN
            WRITE(fp,*) ' --------------------',i2
            WRITE(fp,*) ' OMAP2,4=',obj(iobj)%omap(2),obj(iobj)%omap(4)
            WRITE(fp,*) ' X,Y1   =',x1,y1
            WRITE(fp,*) ' X,Y2   =',x2,y2
          ENDIF
c         Search the wall for intersections:
          s12max = 1.0D+10
          DO iwall = 1, nwall
            x3 = rwall(iwall,1)
            y3 = zwall(iwall,1)
            x4 = rwall(iwall,2)
            y4 = zwall(iwall,2)
            CALL CalcInter(x1,y1,x2,y2,x3,y3,x4,y4,s12,s34) 
            IF (debug) THEN
              WRITE(fp,*) '  CALCINTER :-',i1,i2,iwall
              WRITE(fp,*) '    S12,34  :',s12,s34
              WRITE(fp,*) '    X3,Y3   :',x3,y3
              WRITE(fp,*) '    X4,Y4   :',x4,y4
            ENDIF
            IF (s12.GT.0.0D0.AND.s12.LT.1.0D0.AND.
     .          s34.GT.0.0D0.AND.s34.LT.1.0D0.AND.
     .          s12.LT.s12max) THEN
              s12max = s12
              clist(i1,i2) = iwall
              xlist(i1,i2) = store_x2
              ylist(i1,i2) = store_y2
              IF (debug) WRITE(fp,*) '  *** CUT ***',i1,i2,s12,iwall
            ENDIF
          ENDDO
          IF (clist(i1,i2).EQ.0) THEN
c...        Problem with this cut pair, so delete them from the 
c           list (have to complete the wall by hand at the moment):            
            IF (debug) WRITE(fp,*) ' CUT NOT FOUND, DELETING CUT',i1
            DO i3 = i1, nlist-1
              ilist(i3,:) = ilist(i3+1,:)
            ENDDO
            i1 = i1 - 1
            nlist = nlist - 1
            EXIT
          ENDIF
        ENDDO
      ENDDO

      IF (debug) THEN
        DO i1 = 1, nlist
          WRITE(fp,'(A,4I6,2X,4F14.7)') 
     .      'CLIP LIST:',ilist(i1,:),clist(i1,:),xlist(i1,:),ylist(i1,:)
        ENDDO
      ENDIF

c...  Check if a line segment is cut more than once at either end:
      DO i2 = 1, 2       
        swall = 0
        DO i1 = 1, nlist
          IF (swall(clist(i1,i2)).EQ.0) THEN
            swall(clist(i1,i2)) = 1
          ELSE
            CALL ER('ClipWallToGrid','Wall segment cut at the same '//
     .              'end more than once',*99)
          ENDIF
        ENDDO
      ENDDO

      swall = 0
      DO i1 = 1, nlist
c       Overwrite the approriate end vertices in the wall segments:
        IF (mlist(i1).EQ.1) THEN
          rwall(clist(i1,1),1) = SNGL(xlist(i1,1))
          zwall(clist(i1,1),1) = SNGL(ylist(i1,1))
          rwall(clist(i1,2),2) = SNGL(xlist(i1,2))
          zwall(clist(i1,2),2) = SNGL(ylist(i1,2))
        ELSE
          rwall(clist(i1,1),2) = SNGL(xlist(i1,1))
          zwall(clist(i1,1),2) = SNGL(ylist(i1,1))
          rwall(clist(i1,2),1) = SNGL(xlist(i1,2))
          zwall(clist(i1,2),1) = SNGL(ylist(i1,2))
        ENDIF
c       Select segments that aren't "behind" the targets, and so won't
c       be deleted:
        IF (mlist(i1).EQ.1) THEN
          iwall = clist(i1,1)-1           ! Unlikely that clist(,2) = clist(,1) - 1, 
        ELSE                              ! but if it happens there will be trouble...
          iwall = clist(i1,1)+1           
        ENDIF
        DO WHILE(iwall.NE.clist(i1,2))  
          IF (mlist(i1).EQ.1) THEN
            iwall = iwall + 1 
            IF (iwall.EQ.nwall+1) iwall = 1
          ELSE
            iwall = iwall - 1 
            IF (iwall.EQ.0) iwall = nwall
          ENDIF
          swall(iwall) = 1
          IF (debug) WRITE(fp,*) 'SAVING:',i1,iwall
        ENDDO
      ENDDO

c...  Delete wall segments:
      DO iwall = nwall, 1, -1
        IF (swall(iwall).EQ.1) CYCLE
        IF (delete) THEN
          IF (debug) WRITE(fp,*) 'DELETING:',iwall 
          DO i1 = iwall, nwall-1
            rwall(i1,:) = rwall(i1+1,:)
            zwall(i1,:) = zwall(i1+1,:)
            swall(i1  ) = swall(i1+1  )
          ENDDO
          nwall = nwall - 1
        ELSE
c...      Register that the segment should be deleted but don't remove it:
          IF (debug) WRITE(fp,*) 'MARKING:',iwall 
          rwall(iwall,:) = -9.99D0
          zwall(iwall,:) = -9.99D0
        ENDIF
      ENDDO
      

      IF (debug) THEN
        DO iwall = 1, nwall
          WRITE(fp,'(A,2F14.7,2X,2F14.7)') 
     .      'WALL:',rwall(iwall,1),zwall(iwall,1),
     .              rwall(iwall,2),zwall(iwall,2)
        ENDDO
      ENDIF


      RETURN
 99   WRITE(fp,*) 'I1   =',i1
      WRITE(fp,*) 'I2   =',i2
      WRITE(fp,*) 'S12  =',s12
      WRITE(fp,*) 'IWALL=',iwall
      STOP
      END
c
