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
      SUBROUTINE GenerateTubeGroups
      USE mod_geometry
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      LOGICAL CatchTube
      REAL    ParticleSource

      INTEGER iobj,isrf,it,it1,it2,ic,ic1,ic2,cind1,cind2,cind3,cind4,
     .        itube1,ivtx,ivtx1,ivtx2,in,itsep,region,i,
     .        region_tube(0:100,4),
     .        region_list(0:100,6)
      LOGICAL match
      REAL*8  a1,a2,b1,b2,c1,c2,tab,tcd,rxpt(2),zxpt(2),r0,z0,
     .        volsrc(6)


c      INTEGER, ALLOCATABLE :: cell2obj(:)


      CALL SetupCell2Obj(ncell,IND_CELL)            


c      DO ic = 1, ncell
c          WRITE(0,*) 'CELL3OBJ:',ic,cell2obj(ic)
c      ENDDO


      volsrc = 0.0D0
      region = 1  ! *** LEFT OFF ***
      DO i = 1, 1
        it = 4
        volsrc(region) = volsrc(region) + ParticleSource(it)
              STOP 'sdgsd'
      ENDDO


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
      DO it = 1, ntube
        cind1 = tube(it)%cell_index(LO)
        cind2 = tube(it)%cell_index(HI)
        DO ic1 = cind1, cind2-2
          iobj = cell2obj(ic1)
          isrf = obj(iobj)%iside(3)
          in = 2
          IF (isrf.LT.0) in = 1
          ivtx1 = srf(ABS(isrf))%ivtx(in)
          DO ic2 = ic1+2, cind2
            iobj = cell2obj(ic2)
            isrf = obj(iobj)%iside(3)
            in = 2
            IF (isrf.LT.0) in = 1
            ivtx2 = srf(ABS(isrf))%ivtx(in)
            IF (ivtx1.EQ.ivtx2) THEN
              WRITE(0,*) 'ITUBE!',it
              WRITE(0,*) 'cind1,2=',cind1,cind2
              WRITE(0,*) 'ic1,2  =',ic1,ic2
              WRITE(0,*) 'in     =',in
              itube1 = it
              rxpt(1) = vtx(1,ivtx1)   
              zxpt(1) = vtx(2,ivtx1)   
              itsep = it
              GOTO 10
            ENDIF
          ENDDO
        ENDDO
      ENDDO
 10   CONTINUE

      DO it1 = itube1+1, ntube-1
        cind1 = tube(it1)%cell_index(LO)
        cind2 = tube(it1)%cell_index(HI)
        DO ic1 = cind1, cind2
          iobj = cell2obj(ic1)
          isrf = obj(iobj)%iside(3)
          in = 2
          IF (isrf.LT.0) in = 1
          ivtx1 = srf(ABS(isrf))%ivtx(in)
          DO it2 = it1+1, ntube
            cind3 = tube(it2)%cell_index(LO)
            cind4 = tube(it2)%cell_index(HI)
            DO ic2 = cind3, cind4
              iobj = cell2obj(ic2)
              isrf = obj(iobj)%iside(3)
              in = 2
              IF (isrf.LT.0) in = 1
              ivtx2 = srf(ABS(isrf))%ivtx(in)
              IF (ivtx1.EQ.ivtx2) THEN
                WRITE(0,*) 'ITUBE!',it1,it2
                WRITE(0,*) 'cind1,2=',cind1,cind2
                WRITE(0,*) 'ic1,2  =',ic1,ic2
                WRITE(0,*) 'in     =',in
                rxpt(2) = vtx(1,ivtx1)   
                zxpt(2) = vtx(2,ivtx1)   
                GOTO 20
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
 20   CONTINUE

      WRITE(0,*) 'R,Z0:',r0,z0
      WRITE(0,*) 'XP1 :',rxpt(1),zxpt(1)
      WRITE(0,*) 'XP2 :',rxpt(2),zxpt(2)

      region_tube = 0

c...  High field side SOL rings:
      a1 = r0
      a2 = z0
      b1 = 0.001D0
      b2 = z0
      DO it = itsep, ntube 
        IF (CatchTube(a1,a2,b1,b2,it)) THEN
          WRITE(0,*) 'TUBE CAUGHT HFS:',it
          region_tube(0,1) = region_tube(0,1) + 1
          region_tube(region_tube(0,1),1) = it
        ENDIF
      ENDDO
c...  Identify low field side SOL rings:
      a1 = r0
      a2 = z0
      b1 = r0 + 100.0D0
      b2 = z0
      DO it = itsep, ntube 
        IF (CatchTube(a1,a2,b1,b2,it)) THEN
          WRITE(0,*) 'TUBE CAUGHT LFS:',it
          region_tube(0,2) = region_tube(0,2) + 1
          region_tube(region_tube(0,2),2) = it
        ENDIF
      ENDDO
c...  Secondary PFR:
      a1 = rxpt(2)
      a2 = zxpt(2)
      b1 = 10.0D0 * (rxpt(2) - r0) + r0
      b2 = 10.0D0 * (zxpt(2) - z0) + z0
      DO it = itsep, ntube 
        IF (CatchTube(a1,a2,b1,b2,it)) THEN
          WRITE(0,*) 'TUBE CAUGHT SECONDARY PFR:',it
          region_tube(0,3) = region_tube(0,3) + 1
          region_tube(region_tube(0,3),3) = it
        ENDIF
      ENDDO
c...  Primary PFR:
      a1 = rxpt(1)
      a2 = zxpt(1)
      b1 = 10.0D0 * (rxpt(1) - r0) + r0
      b2 = 10.0D0 * (zxpt(1) - z0) + z0
      DO it = itsep, ntube 
        IF (CatchTube(a1,a2,b1,b2,it)) THEN
          WRITE(0,*) 'TUBE CAUGHT PRIMARY   PFR:',it
          region_tube(0,4) = region_tube(0,4) + 1
          region_tube(region_tube(0,4),4) = it
        ENDIF
      ENDDO

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
              STOP 'sdgsd'
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
      SUBROUTINE ProcessGrid
      USE mod_geometry
      USE mod_legacy
      IMPLICIT none

      INTEGER status, grd_load_method

      INTEGER, PARAMETER :: GRD_LOAD_NEW = 2, GRD_LOAD_OLD = 1

      grd_load_method = 2


      SELECTCASE (grd_load_method)
        CASE (GRD_LOAD_NEW)
          CALL LoadGeneralisedGrid
        CASE (GRD_LOAD_OLD)
          CALL LoadObjects('osm_geometry.raw',status)  ! Change to LoadGeometryObjects...
          IF (status.EQ.-1) 
     .      CALL ER('SetupGrid','Geometry data not found',*99)

          CALL LoadGrid('osm.raw')

          CALL LoadLegacyData('osm_legacy.raw')
        CASE DEFAULT
          CALL ER('ProcessGrid','Unrecognized grid source option',*99)
      ENDSELECT
 
      RETURN
 99   STOP
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
      RECURSIVE SUBROUTINE FindGridCell(nknot,knot,NUMZONE,izone,
     .                                   condition,index1,index2)
      USE mod_grid
      IMPLICIT none

      INTEGER nknot,index1,index2,NUMZONE,izone(NUMZONE+1,NUMZONE),
     .        condition
      TYPE(type_grid_cell) :: knot(0:nknot)      

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
            IF     (condition.EQ.1.AND.
     .              ABS(knot(i1)%rv(1)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(1)-knot(i2)%zv(1)).LT.DTOL) THEN

              if (output) then
                 WRITE(0,*) 'XPOINT SOL:',i1,i2
                 WRITE(6,*) 'XPOINT SOL:',i1,i2
              endif
                 
              index2 = i2
              RETURN
            ELSEIF     (condition.EQ.6.AND.
     .              ABS(knot(i1)%rv(2)-knot(i2)%rv(2)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(2)-knot(i2)%zv(2)).LT.DTOL) THEN
c... *** WHAT'S THE POINT OF THIS? ***
              STOP 'MYSTERIOUS condition.EQ.6'
              if (output) then
                 WRITE(0,*) 'XPOINT CORE/PFZ:',i1,i2
                 WRITE(6,*) 'XPOINT CORE/PFZ:',i1,i2
              endif
                 
              index2 = i2
              RETURN

            ELSEIF (condition.EQ.2.AND.
     .              ABS(knot(i1)%rv(1)-knot(i2)%rv(2)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(1)-knot(i2)%zv(2)).LT.DTOL.AND.
     .              ABS(knot(i1)%rv(4)-knot(i2)%rv(3)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(4)-knot(i2)%zv(3)).LT.DTOL) THEN

              if (output) then 
                WRITE(0,*) 'SIDE INWARD 41:',i1,i2
                WRITE(6,*) 'SIDE INWARD 41:',i1,i2
              endif

              index2 = i2
              RETURN

            ELSEIF (condition.EQ.3.AND.
     .              ABS(knot(i1)%rv(3)-knot(i2)%rv(2)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(3)-knot(i2)%zv(2)).LT.DTOL.AND.
     .              ABS(knot(i1)%rv(4)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(4)-knot(i2)%zv(1)).LT.DTOL) THEN
              IF (output) THEN
                WRITE(0,*) 'SIDE UP 34:',i1,i2
                WRITE(6,*) 'SIDE UP 34:',i1,i2
              ENDIF
              index2 = i2
              RETURN
            ELSEIF (condition.EQ.4.AND.
     .              ABS(knot(i1)%rv(2)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(2)-knot(i2)%zv(1)).LT.DTOL.AND.
     .              ABS(knot(i1)%rv(3)-knot(i2)%rv(4)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(3)-knot(i2)%zv(4)).LT.DTOL) THEN


               if (output) then 
                 WRITE(0,*) 'SIDE OUTWARD 23:',i1,i2
                 WRITE(6,*) 'SIDE OUTWARD 23:',i1,i2
               endif

              index2 = i2
              RETURN

            ELSEIF (condition.EQ.5.AND.
     .              ABS(knot(i1)%rv(1)-knot(i2)%rv(4)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(1)-knot(i2)%zv(4)).LT.DTOL.AND.
     .              ABS(knot(i1)%rv(2)-knot(i2)%rv(3)).LT.DTOL.AND.
     .              ABS(knot(i1)%zv(2)-knot(i2)%zv(3)).LT.DTOL) THEN
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
      USE mod_grid
      IMPLICIT none

c...  Input:
      INTEGER grdfp

c...  Output:
      INTEGER irsep,irsep2,irwall,irtrap,nrs,ikti,ikto,nks(100)
      INTEGER ik,ir
      LOGICAL connected
      
      INTEGER, PARAMETER :: NUMZONE = 5
      REAL*8,  PARAMETER :: HI   = 1.0E+20
      REAL*8,  PARAMETER :: DTOL = 1.0D-07

      INTEGER   nknot,i1,i2,z1,r1,kind,nxpt,ixpt(0:2,2),cxpt(0:2,2),i3,
     .          i4,izone(NUMZONE+1,NUMZONE),newi1,icore(0:2,2),id,
     .          tmpnks,ir1,istart,fp,maxik,maxir,outfp,count

      LOGICAL   cont,deleteknot,output,swap,cell_deletion
      REAL*8    vrmin,vzmin,vrmax,vzmax,rspan,zspan,area
      REAL*8    rvdp(4),zvdp(4),areadp
      CHARACTER buffer*1000

      INTEGER, ALLOCATABLE :: imap(:,:)

c     TYPE type_cell
c     INTEGER :: index,ik,ir,nv,rzone,zzone,xpt,map
c     REAL    :: rcen,zcen,bratio,rv(4),zv(4)
c     ENDTYPE type_cell
      
      TYPE(type_grid_cell),ALLOCATABLE :: knot(:)
c     
c     jdemod - flag for turning on and off debugging output
c     
      output = .TRUE.

      outfp = 0

c...  Read the knot data:

      grd_format = 1
      grd_filename = 'iterm.carre.105'
      grdfp = 99
      OPEN(UNIT=grdfp,FILE=TRIM(grd_filename),ACCESS='SEQUENTIAL',
     .     ERR=98)     

      SELECTCASE (grd_format)
c     ----------------------------------------------------------------
      CASE (GRD_FORMAT_SONNET)
c...    Find the start of the cell/knot information in the grid file:
        WRITE(buffer,'(1000X)')
        DO WHILE (buffer(4:8).NE.'=====')
           READ(grdfp,'(A10)',END=98) buffer
        ENDDO
c...    Scan the file to see how many cells are in the grid:
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
 10     CONTINUE  ! EOF
        nknot = nknot + 1
c...    Setup arrays:
        ALLOCATE(knot(0:nknot))
        ALLOCATE(imap(maxik,0:maxir))
        nknot = 0
c...    Load grid:
        REWIND(grdfp)
        WRITE(buffer,'(1000X)')
        DO WHILE (buffer(4:8).NE.'=====')
           READ(grdfp,'(A10)',END=98) buffer
        ENDDO
        DO WHILE(nknot.EQ.0.OR.buffer(4:10).EQ.'Element')
c          READ(grdfp,'(A50)',END=10) buffer
c          WRITE(0,*) 'BUFFER:',buffer(1:50)
c          BACKSPACE(grdfp)
          nknot = nknot + 1
          READ(grdfp,80,END=97) knot(nknot)%index,
     .                          knot(nknot)%ik   ,knot(nknot)%ir, 
     .                          knot(nknot)%rv(2),knot(nknot)%zv(2),
     .                          knot(nknot)%rv(3),knot(nknot)%zv(3)
c          WRITE(0,*) 'VAL:',knot(nknot)%rv(2)
          READ(grdfp,81,END=97) knot(nknot)%bratio,
     .                          knot(nknot)%rcen ,knot(nknot)%zcen
          READ(grdfp,82,END=97) knot(nknot)%rv(1),knot(nknot)%zv(1),
     .                          knot(nknot)%rv(4),knot(nknot)%zv(4)
          knot(nknot)%nv = 4
c...      Dividing line in grid file:       
          READ(grdfp,*)
          READ(grdfp,'(A10)',END=20) buffer
          BACKSPACE(grdfp)
        ENDDO
 80     FORMAT(10X,I5,4X,I3,1x,I3,4x,E17.10,1X,E17.10,8X,E17.10,
     .         1X,E17.10)
 81     FORMAT(18X,E17.10,14x,E17.10,1x,E17.10)
 82     FORMAT(30X,E17.10,1X,E17.10,8x,E17.10,1X,E17.10)
 20     CONTINUE  ! EOF
c     ----------------------------------------------------------------
      CASE DEFAULT
        CALL ER('LoadGeneralisedGrid','Unknown grid type',*99)
      ENDSELECT

      IF (output) THEN
        WRITE(outfp,*) 'GRID LOADED'
        WRITE(outfp,*) '  NKNOT=',nknot
        WRITE(outfp,*) '  MAXIK=',maxik
        WRITE(outfp,*) '  MAXIR=',maxir
      ENDIF
 
c...  Delete zero volume cells:
c     
      IF (output) THEN
        WRITE(outfp,*) 'REMOVING ZERO VOLUME CELLS'
      ENDIF
      cell_deletion = .FALSE.
      i1 = 1
      DO WHILE(i1.LE.nknot)
        IF     (.FALSE..AND.
     .       (knot(i1)%ik.EQ.0.OR.knot(i1)%ik.EQ.maxik-1.OR. ! virtual cells on end of rings
     .        knot(i1)%ir.EQ.0.OR.knot(i1)%ir.EQ.maxir-1)) THEN ! virtual rings
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
          IF (output) THEN
            WRITE(outfp,*) ' >',i1,knot(i1)%ik,knot(i1)%ir
          ENDIF
          IF (i1.LT.nknot) THEN
            DO i2 = i1, nknot
              CALL MoveGridCell(knot(i2+1),knot(i2))
            ENDDO
          ENDIF
          nknot = nknot - 1
        ENDIF
      ENDDO
      IF (output) THEN
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
c       Not implemented yet

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
      IF (output) THEN
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

      IF (output) THEN
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
        CALL FindGridCell(nknot,knot,NUMZONE,izone,1,i1,i2)
        IF (i2.NE.-1) THEN
          nxpt = nxpt + 1
          ixpt(nxpt,1) = i1
          ixpt(nxpt,2) = i2
          knot(i1)%xpt = i2
          knot(i2)%xpt = i1
        ENDIF
      ENDDO
c...  No x-points found, for some reason:
      IF (output) THEN
        WRITE(outfp,*) 'X-POINT CELLS:'
        WRITE(outfp,*) '  IXPT1=',ixpt(1:nxpt,1)
        WRITE(outfp,*) '  IXPT2=',ixpt(1:nxpt,2)
      ENDIF
      IF (nxpt.EQ.0) 
     .  CALL ER('LoadGeneralisedgrid','No x-points found',*99)
      IF (nxpt.GT.2) 
     .  CALL ER('LoadGeneralizedGrid','More than 2 x-points '//
     .          'identified, not good',*99)

c...  
      cxpt = 0
      icore = 0
      DO i4 = 1, nxpt
        DO i3 = 1, 2
          newi1 = ixpt(i4,i3)
          cont = .TRUE.
          DO WHILE (cont)
            i1 = newi1
            cxpt(i4,i3) = cxpt(i4,i3) + 1
            cont = .FALSE.
            CALL FindGridCell(nknot,knot,NUMZONE,izone,2,i1,i2)
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
            CALL FindGridCell(nknot,knot,NUMZONE,izone,3,i1,i2)  ! Move forward on a ring
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
     .               'processing x-points -- malformed grid?',*99)
        ENDDO
c...    Swap the cells, if necessary:
        IF (icore(i4,1).LT.0) THEN
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
      IF (output) THEN
        WRITE(outfp,*) 'X-POINT SCANNING:' 
        WRITE(outfp,*) '   IXTP1 =',ixpt (1:nxpt,1)
        WRITE(outfp,*) '   CXPT1 =',cxpt (1:nxpt,1)
        WRITE(outfp,*) '   ICORE1=',icore(1:nxpt,1)
        WRITE(outfp,*) 
        WRITE(outfp,*) '   IXTP2 =',ixpt (1:nxpt,2)
        WRITE(outfp,*) '   CXPT2 =',cxpt (1:nxpt,2)
        WRITE(outfp,*) '   ICORE2=',icore(1:nxpt,2)
      ENDIF

c...  Check that the x-points are ordered properly, with the primary x-point
c     at index 1, and whether or not the double-null grid is connected:
      connected = .FALSE.
      IF (nxpt.GT.1) THEN
        swap = .FALSE.
        IF     (nxpt.GT.1.AND.cxpt(1,1).EQ.cxpt(2,1)) THEN
c...      Connected grid indicated:
          IF (knot(ixpt(1,1))%zcen.GT.0.0) swap = .TRUE.  ! Want lower x-point to be "primary"
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
        ENDIF
      ENDIF

c...  Assemble the grid:
c     ------------------------------------------------------------------

c...  Location of the primary separatrix is known:
      irsep  = cxpt(1,1)
      irsep2 = irsep

c...  Build the core rings:
c     ------------------------------------------------------------------
      IF (output) WRITE(outfp,*) 'PROCESSING CORE RINGS'

      nrs = 1
      ik = 1
      ir = 1
      i1 = icore(1,1) 
      DO WHILE(ir.LT.irsep)
        imap(1,ir) = i1
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
          CALL FindGridCell(nknot,knot,NUMZONE,izone,3,i1,i2)  ! Move forward
          IF (i2.NE.-1) THEN
            IF (i2.NE.imap(1,ir)) THEN
              i1 = i2 
              ik = ik + 1
              imap(ik,ir) = i1
              cont = .TRUE.
c             IF (output) WRITE(0,*) 'CORE MAP:',ik,ir,i1
            ELSE
c             Ring has closed on itself, finished here
            ENDIF
          ELSE
            CALL ER('LoadGeneralisedgrid','Cell sequence broken in '//
     .              'the core',*99)
          ENDIF
        ENDDO
        nks(ir) = ik
c...    Step outward to the next ring, still in the core:        
        CALL FindGridCell(nknot,knot,NUMZONE,izone,4,imap(1,ir),i2)  ! And to the left
        IF (i2.NE.-1) THEN        
          i1 = i2
          ik = 1
          ir = ir + 1
        ELSE
          CALL ER('LoadGeneralisedgrid','Radial step has failed',*99)
        ENDIF
      ENDDO
      nrs = ir
      IF (output) THEN
        WRITE(outfp,*) ' CORE RINGS:'
        WRITE(outfp,*) '  NRS  =',nrs
        WRITE(outfp,*) '  IRSEP=',irsep
        WRITE(outfp,*) '  NKS  =',nks(1:irsep-1)
      ENDIF

c...  SOL rings:
c     ------------------------------------------------------------------
      IF (output) WRITE(outfp,*) 'PROCESSING SOL RINGS'

c...  Step out of the core:
      i1 = imap(1,irsep-1)
      CALL FindGridCell(nknot,knot,NUMZONE,izone,4,i1,i2)
      IF (i2.NE.-1) THEN  
c...    We should be back at the primary x-point:
        IF (i2.NE.ixpt(1,1)) 
     .    CALL ER('LoadGeneralisedGrid','Not back at the x-point '//
     .            'cell for some reason',*99)
c...    Move backward along the ring to find the target cell:
        i1 = i2
        DO WHILE(i2.NE.-1)
          CALL FindGridCell(nknot,knot,NUMZONE,izone,5,i1,i2) ! Moving backward 
          IF (i2.NE.-1) i1 = i2
        ENDDO
      ELSE
        CALL ER('LoadGeneralisedgrid','Radial step out from the '//
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
        CALL FindGridCell(nknot,knot,NUMZONE,izone,3,i1,i2)
        IF (i2.NE.-1) THEN
          i1 = i2 
          ik = ik + 1
          imap(ik,ir) = i1
          cont = .TRUE.
          IF (output) WRITE(0,*) 'SOL MAP (INNER):',ik,ir,i1
        ENDIF
        IF (.NOT.cont) THEN
c...      Step radially outward if ring is finished, i.e. a cell has been reached
c         for which no "next" cell is found when moving along the ring:
          nks(ir) = ik
          i1 = imap(1,ir)
          CALL FindGridCell(nknot,knot,NUMZONE,izone,4,i1,i2)  ! Step radially          
          IF (i2.NE.-1) THEN
             i1 = i2
             ik = 1
             ir = ir + 1
             imap(ik,ir) = i1
             cont = .TRUE.
          ELSE
            CALL ER('LoadGeneralisedgrid','Radial step '//
     .              'failed in the SOL',*99)
          ENDIF
        ENDIF
      ENDDO
      irwall = ir
      irtrap = ir + 1
      nrs = ir

c...  SOL rings associated with double-null grids:
c     ------------------------------------------------------------------
      IF (nxpt.GT.1) THEN
c...    Find scrape-off layer rings that were not processed yet, which 
c       happens with double-null grids:
        IF (output) WRITE(outfp,*) 'DOUBLE NULL GRID, PROCESSING '// 
     .                             'REMAINDER OF SOL'
c...    Register all cells that have been mapped to the grid:
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            knot(imap(ik,ir))%map = 1
          ENDDO
        ENDDO

c...    Start from the last cell on the last core ring:
        i1=imap(nks(irsep-1),irsep-1) 
c        i1 = imap(1,irsep-1)
c        CALL FindGridCell(nknot,knot,NUMZONE,izone,5,i1,i2)
c        IF (i2.EQ.-1) 
c     .       CALL ER('LoadGeneralisedgrid','Core map problems',*99)
c        i1 = i2

c...    Move radially outward into the SOL:
        CALL FindGridCell(nknot,knot,NUMZONE,izone,4,i1,i2)
        IF (i2.NE.-1) THEN
          i1 = i2
          IF (knot(i2)%map.EQ.1) THEN  
c...        Keep going until a cell with no assigned mapping is found:
            cont = .TRUE.
            DO WHILE(cont)
              cont = .FALSE.
              CALL FindGridCell(nknot,knot,NUMZONE,izone,4,i1,i2)
              IF (i2.NE.-1) THEN 
                i1 = i2
                IF (knot(i1)%map.NE.0) cont = .TRUE.
              ELSE

              ENDIF
            ENDDO
          ENDIF
        ELSE
           CALL ER('LoadGeneralisedGrid','Bad IR step to SOL',*99)
        ENDIF


c *** LEFT OFF ***

c...  An unmapped cell has been found, proceed to target:
         cont = .TRUE.
         DO WHILE(cont)
            cont = .FALSE.
            CALL FindGridCell(nknot,knot,NUMZONE,izone,5,i1,i2)
            IF (i2.NE.-1) THEN 
               i1 = i2
               cont = .TRUE.
            ENDIF
         ENDDO

c...  Target located, start mapping the low field SOL:
         ik = 1
         ir = nrs + 1
         IF (connected) irsep2 = ir
         imap(ik,ir) = i1
         cont = .TRUE.
         DO WHILE(cont)
            cont = .FALSE.
c...  Move along the ring:
            CALL FindGridCell(nknot,knot,NUMZONE,izone,3,i1,i2)
            IF (i2.NE.-1) THEN
               i1 = i2 
               ik = ik + 1
               imap(ik,ir) = i1
               cont = .TRUE.
c     
c     jdemod
               IF (output) then 
                  WRITE(0,'(a,4i6,5(1x,g12.5))') 'OUTER SOL MAP:',
     >                 ik,ir,i1,imap(1,ir),
     >                 knot(i1)%rcen,knot(i1)%zcen
                  WRITE(6,'(a,4i6,5(1x,g12.5))') 'OUTER SOL MAP:',
     >                 ik,ir,i1,imap(1,ir),
     >                 knot(i1)%rcen,knot(i1)%zcen
               endif
c     IF (output) WRITE(0,*) 'OUTER SOL MAP:',ik,ir,i1
            ENDIF
c...  Step radially outward if ring is finished:
            IF (.NOT.cont) THEN
c     
c     jdemod
               IF (output) then 
                  WRITE(0,*) 'STEPPING OUT:',ik,ir,i1
                  WRITE(6,*) 'STEPPING OUT:',ik,ir,i1
               endif
c     IF (output) WRITE(0,*) 'STEPPING OUT:',ik,ir,i1
               nks(ir) = ik
               i1 = imap(1,ir)
               CALL FindGridCell(nknot,knot,NUMZONE,izone,4,i1,i2)          
               IF (i2.NE.-1) THEN
                  i1 = i2
                  ik = 1
                  ir = ir + 1
                  imap(ik,ir) = i1
                  cont = .TRUE.
c     
c     jdemod
                  IF (output) then 
                     WRITE(0,*) 'OUTER SOL MAP NEW RING:',ik,ir,i1
                     WRITE(6,*) 'OUTER SOL MAP NEW RING:',ik,ir,i1
                  endif
c     IF (output) WRITE(0,*) 'OUTER SOL MAP NEW RING:',ik,ir,i1
               ELSE
c...  Assume the outer boundary of the grid:
c     
c     jdemod
                  IF (output) then 
                     WRITE(0,*) 'ASSUMING OUTER GRID BOUNDARY'
                     WRITE(6,*) 'ASSUMING OUTER GRID BOUNDARY'
                  endif
c     IF (output) WRITE(0,*) 'ASSUMING OUTER GRID BOUNDARY'
               ENDIF
            ENDIF
         ENDDO
         irwall = ir
         irtrap = ir + 1
         nrs = ir

c...  Register all knots that have been mapped to the grid:
         DO ir = 1, nrs
            DO ik = 1, nks(ir)
               knot(imap(ik,ir))%map = 1
            ENDDO
         ENDDO

c     
c     jdemod
         IF (output) then 
            WRITE(0,*) 'PROCESSING SECONDARY PFZ'
            WRITE(6,*) 'PROCESSING SECONDARY PFZ'
         endif
c     IF (output) WRITE(0,*) 'PROCESSING SECONDARY PFZ'
c...  Process the secondary x-point PFR, which is just considered part of the
c     SOL for generalized grids:
         i1 = knot(ixpt(2,1))%xpt
c     i1 = ixpt(2)
c...  Move into the PFR:
         CALL FindGridCell(nknot,knot,NUMZONE,izone,2,i1,i2)
         IF (i2.EQ.-1) 
     .        CALL ER('LoadGeneralisedgrid','PFR2 problems',*99)
         i1 = i2
         ik = 1
         ir = nrs + 1
         imap(ik,ir) = i1

c...  Proceed to target:
c     
c     jdemod
         IF (output) then 
            WRITE(0,*) 'LOOKING FOR TARGET'
            WRITE(6,*) 'LOOKING FOR TARGET'
         endif
c     IF (output) WRITE(0,*) 'LOOKING FOR TARGET'
         cont = .TRUE.
         DO WHILE(cont)
            cont = .FALSE.
            CALL FindGridCell(nknot,knot,NUMZONE,izone,5,i1,i2)
c     
c     jdemod
            IF (output) then 
               WRITE(0,*) 'MOVING',i1,i2,istart
               WRITE(6,*) 'MOVING',i1,i2,istart
            endif
c     IF (output) WRITE(0,*) 'MOVING',i1,i2,istart
            IF (i2.NE.-1) THEN 
               i1 = i2
               ik = ik + 1
               imap(ik,ir) = i1
               cont = .TRUE.
            ENDIF
         ENDDO
         
c...  Target located, start mapping the secondary PFR:
c     
c     jdemod
         IF (output) then 
            WRITE(0,*) 'TARGET LOCATED'
            WRITE(6,*) 'TARGET LOCATED'
         endif
c     IF (output) WRITE(0,*) 'TARGET LOCATED'
         ik = 1
         ir = nrs + 1
         imap(ik,ir) = i1
         cont = .TRUE.
         DO WHILE(cont)
            cont = .FALSE.
c...  Move along the ring:
            CALL FindGridCell(nknot,knot,NUMZONE,izone,3,i1,i2)
c     
c     jdemod
            IF (output) then 
               WRITE(0,*) 'MOVING',i1,i2
               WRITE(6,*) 'MOVING',i1,i2
            endif
c     IF (output) WRITE(0,*) 'MOVING',i1,i2
            IF (i2.NE.-1) THEN
               i1 = i2 
               ik = ik + 1
               imap(ik,ir) = i1
               cont = .TRUE.
c     
c     jdemod
               IF (output) then 
                  WRITE(0,*) '2ND PFR MAP:',ik,ir,i1,knot(i1)%zcen
                  WRITE(6,*) '2ND PFR MAP:',ik,ir,i1,knot(i1)%zcen
               endif
c     IF (output) WRITE(0,*) '2ND PFR MAP:',ik,ir,i1,knot(i1)%zcen
            ENDIF
c...  Step radially outward if ring is finished:
            IF (.NOT.cont) THEN
               nks(ir) = ik
               i1 = imap(1,ir)
               CALL FindGridCell(nknot,knot,NUMZONE,izone,2,i1,i2)          
               IF (i2.NE.-1) THEN
                  i1 = i2
                  ik = 1
                  ir = ir + 1
                  imap(ik,ir) = i1
                  cont = .TRUE.
c     
c     jdemod
                  IF (output) then  
                     WRITE(0,*) '2ND PFR MAP NEW RING:',ik,ir,i1,
     >                    knot(i1)%zcen
                     WRITE(6,*) '2ND PFR MAP NEW RING:',ik,ir,i1,
     >                    knot(i1)%zcen
                  endif
c     IF (output) 
c     .          WRITE(0,*) '2ND PFR MAP NEW RING:',ik,ir,i1,
c     .          knot(i1)%zcen
               ELSE
c...  Assume the outer boundary of the grid:
               ENDIF
            ENDIF
         ENDDO
         irwall = ir
         irtrap = ir + 1
         nrs = ir

      ENDIF                     ! Done processing double-null rings


c...  Primary private flux zone (PFZ) rings:
c     ------------------------------------------------------------------
c     
      IF (output) WRITE(outfp,*) 'PROCESSING PRIMARY PFZ'


c...  Process the primary x-point PFR:
      i1 = knot(ixpt(1,1))%xpt
c...  Move into the PFR:
      CALL FindGridCell(nknot,knot,NUMZONE,izone,2,i1,i2)
      IF (i2.EQ.-1) 
     .     CALL ER('LoadGeneralisedgrid','PFR1 problems',*99)
c     
c     jdemod
      IF (output) then
         WRITE(0,*) '  INTO PFZ',knot(i1)%index,knot(i2)%index
         WRITE(6,*) '  INTO PFZ',knot(i1)%index,knot(i2)%index
      endif
c     IF (output) 
c     .  WRITE(0,*) '  INTO PFZ',knot(i1)%index,knot(i2)%index
      i1 = i2
c...  Proceed to target:
      cont = .TRUE.
      DO WHILE(cont)
         cont = .FALSE.
         CALL FindGridCell(nknot,knot,NUMZONE,izone,5,i1,i2)
c     jdemod
c     IF (output) then 
c     WRITE(0,*) '  TO TARGET',i1,i2
c     WRITE(0,*) '  TO TARGET',knot(i1)%index
c     WRITE(0,*) '  TO TARGET',knot(i2)%index
c     WRITE(6,*) '  TO TARGET',i1,i2
c     WRITE(6,*) '  TO TARGET',knot(i1)%index
c     WRITE(6,*) '  TO TARGET',knot(i2)%index
c     endif
c     IF (output) WRITE(0,*) '  TO TARGET',i1,i2
c     IF (output) WRITE(0,*) '  TO TARGET',knot(i1)%index
c     IF (output) WRITE(0,*) '  TO TARGET',knot(i2)%index
c     STOP 'tet'
         IF (i2.NE.-1) THEN 
            i1 = i2
            cont = .TRUE.
         ENDIF
      ENDDO
c...  Target located, start mapping the primary PFR:
c     
c     jdemod
      IF (output) then 
         WRITE(0,*) 'PROCESSING PRIMARY PFZ, TARGET LOCATED'
         WRITE(6,*) 'PROCESSING PRIMARY PFZ, TARGET LOCATED'
      endif
c     IF (output) WRITE(0,*) 'PROCESSING PRIMARY PFZ, TARGET LOCATED'
      ik = 1
      ir = nrs + 1
      imap(ik,ir) = i1
      cont = .TRUE.
      DO WHILE(cont)
         cont = .FALSE.
c...  Move along the ring:
         CALL FindGridCell(nknot,knot,NUMZONE,izone,3,i1,i2)
         IF (i2.NE.-1) THEN
            i1 = i2 
            ik = ik + 1
            imap(ik,ir) = i1
            cont = .TRUE.
c     jdemod
c     IF (output) then 
c     WRITE(0,*) '1ST PFR MAP:',ik,ir,i1,knot(i1)%zcen
c     WRITE(6,*) '1ST PFR MAP:',ik,ir,i1,knot(i1)%zcen
c     endif
c     IF (output) WRITE(0,*) '1ST PFR MAP:',ik,ir,i1,knot(i1)%zcen
         ENDIF
c...  Step radially outward if ring is finished:
         IF (.NOT.cont) THEN
            nks(ir) = ik
            i1 = imap(1,ir)
            CALL FindGridCell(nknot,knot,NUMZONE,izone,2,i1,i2)          
            IF (i2.NE.-1) THEN
               i1 = i2
               ik = 1
               ir = ir + 1
               imap(ik,ir) = i1
               cont = .TRUE.
c     
c     jdemod
               IF (output) then 
                  WRITE(0,'(a,4i6,5(1x,g12.5))') 'PRIM PFZ MAP:',
     >                 ik,ir,i1,imap(1,ir),
     >                 knot(i1)%rcen,knot(i1)%zcen
                  WRITE(6,'(a,4i6,5(1x,g12.5))') 'PRIM PFZ MAP:',
     >                 ik,ir,i1,imap(1,ir),
     >                 knot(i1)%rcen,knot(i1)%zcen
               endif
c     IF (output) 
c     .        WRITE(0,*) '1ST PFR MAP NEW RING:',ik,ir,i1,knot(i1)%zcen
            ELSE
c...  Assume the outer boundary of the grid:
            ENDIF
c     jdemod
            IF (output) then 
               WRITE(0,*) 'PROCESSING PRIMARY PFZ, BUZZING...'
               WRITE(6,*) 'PROCESSING PRIMARY PFZ, BUZZING...'
            endif
c     IF (output) WRITE(0,*) 'PROCESSING PRIMARY PFZ, BUZZING...'
         ENDIF
      ENDDO
      nrs = ir

c...  Need to reorder the rings in the primary PFZ:
      DO i1 = 0, (nrs-irtrap+1)/2-1
c     
c     jdemod
         IF (output) then 
            WRITE(0,*) 'I1???=',i1
            WRITE(6,*) 'I1???=',i1
         endif
c     IF (output) WRITE(0,*) 'I1???=',i1
         tmpnks = nks(irtrap+i1)
         DO ik = 1, nks(irtrap+i1)
            imap(ik,0) = imap(ik,irtrap+i1)
         ENDDO
         nks(irtrap+i1) = nks(nrs-i1)
         DO ik = 1, nks(nrs-i1)
            imap(ik,irtrap+i1) = imap(ik,nrs-i1)
         ENDDO
         nks(nrs-i1) = tmpnks
         DO ik = 1, tmpnks
            imap(ik,nrs-i1) = imap(ik,0)
         ENDDO
      ENDDO

c...  Find IKTO,IKTI:
      ikto = -1
      ikti = -1
      DO ik = 1, nks(irsep)
         IF (connected) THEN
            IF (imap(ik,irsep).EQ.ixpt(1,1)          ) ikto = ik ! Not sure this will always work...
            IF (imap(ik,irsep).EQ.knot(ixpt(2,1))%xpt) ikti = ik - 1
c     IF (imap(ik,irsep).EQ.ixpt(1,1)) ikto = ik      ! Not sure this will always work...
c     IF (imap(ik,irsep).EQ.ixpt(2,1)) ikti = ik - 1
c     IF (imap(ik,irsep).EQ.knot(ixpt(2,1))%xpt) ikti = ik - 1  ! Not sure this will always work...
         ELSE
            IF (imap(ik,irsep).EQ.ixpt(1,1)          ) ikto = ik - 1
            IF (imap(ik,irsep).EQ.knot(ixpt(1,1))%xpt) ikti = ik 
         ENDIF
      ENDDO
      IF (ikto.EQ.-1.OR.ikti.EQ.-1)
     .     CALL ER('LoadGeneralisedgrid','IKTI or IKTO not found',*99)


      RETURN
 96   CALL ER('LoadGeneralisedgrid','Grid file not found',*99)
 97   CALL ER('LoadGeneralisedgrid','Unexpected end-of-file',*99)
 98   CALL ER('LoadGeneralisedgrid','Problem accessing grid file',*99)
 99   WRITE(0,*) 'GRID:',TRIM(grd_filename)
      WRITE(0,*) 'IXPT:',ixpt(1,1),ixpt(2,1)
      STOP  
      END
c
