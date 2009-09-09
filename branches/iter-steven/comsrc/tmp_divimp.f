
c
c ======================================================================
c
c subroutine: CalcTubeDimentions
c
      SUBROUTINE CalcTubeDimensions(tube_3D_data,dangle)
      IMPLICIT none

c      SUBROUTINE CalcTubeDimensions(xin,yin,zin,
c     .             MAXSURFACE,MAXPOINTS,nsur,npts,hsur,vsur)
c      IMPLICIT none
c
c      INTEGER, INTENT(IN) :: MAXSURFACE,MAXPOINTS
c      REAL   , INTENT(IN) :: xin,yin,zin
c      INTEGER, INTENT(OUT) :: nsur,npts(MAXSURFACE),hsur(MAXSURFACE)
c      REAL*8 , INTENT(OUT) :: vsur(3,MAXPOINTS,0:MAXSURFACE)

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      REAL tube_3D_data(5,MAXNKS,MAXNRS),dangle

      INTEGER ir,isegment,nsegment,mode,ncell,icell
      REAL    rho_start,rho_end,rho1,rho2,length_total

      TYPE    type_tube_3D
         REAL      :: dangle
         INTEGER   :: ncell
         REAL*8    :: area  (1000)
         REAL*8    :: length(1000)
         REAL*8    :: volume(1000)
         INTEGER   :: n
         INTEGER   :: index(  200000)
         REAL*8    :: v    (3,200000)
      ENDTYPE type_tube_3D

      TYPE(type_tube_3D) :: tube_3D(0:5)

      WRITE(0,*) '---------------------------------'
      WRITE(0,*) 'NOT CALCULATING TUBE DIMENTSIONS!'
      WRITE(0,*) '---------------------------------'
      RETURN


      nsegment = 25

      DO ir = 5, 9 !  irsep, irsep+1

        rho_start = rho(ir,IN14 ) * 1.0001
        rho_end   = rho(ir,OUT23) * 0.9999

        mode = 1

        DO isegment = 1, nsegment

          rho1 = rho_start + REAL(isegment-1) / REAL(nsegment) * 
     .                       (rho_end - rho_start)
          rho2 = rho_start + REAL(isegment  ) / REAL(nsegment) * 
     .                       (rho_end - rho_start)

          CALL ProcessFluxTube(rho1,rho2,REAL(ir),mode,tube_3D(1:5))
c          CALL ProcessFluxTube(rho1,rho2,REAL(ir),
c     .           mode,tube_3D(1:5),
c     .           MAXSURFACE,MAXPOINTS,nsur,npts,hsur,vsur)

          IF (mode.EQ.1) THEN
            tube_3D(0) = tube_3D(1)
            ncell = tube_3D(0)%ncell
          ELSE
            tube_3D(0)%area   = tube_3D(0)%area   + tube_3D(1)%area   
            tube_3D(0)%length = tube_3D(0)%length + tube_3D(1)%length 
            tube_3D(0)%volume = tube_3D(0)%volume + tube_3D(1)%volume
          ENDIF

          mode = 2
        ENDDO

        tube_3D(0)%length = tube_3D(0)%length / REAL(nsegment) 

        dangle = tube_3D(0)%dangle

        length_total = 0.0
        DO icell = 1, ncell

          tube_3D_data(1,icell,ir) = length_total

          length_total = length_total + SNGL(tube_3D(0)%length(icell))

          WRITE(6,'(A,2I6,2F10.5,2X,2F10.5)') 
     .      'TOTAL VOLUME:',icell,ir,
     .      SNGL(tube_3D(0)%volume(icell) *
     .      (360.0D0/tube_3D(0)%dangle)),
     .      kvols(icell,ir),
     .      length_total,
     .      ksb(icell,ir)

          tube_3D_data(2,icell,ir) = length_total
          tube_3D_data(3,icell,ir) = length_total - 
     .                    0.5 * SNGL(tube_3D(0)%length(icell))       
          tube_3D_data(4,icell,ir) = tube_3D(0)%area  (icell)
          tube_3D_data(5,icell,ir) = tube_3D(0)%volume(icell)

        ENDDO
      ENDDO

      RETURN
 99   STOP
      END
c     
c ======================================================================
c
c subroutine: SelectGridRegoin_DIVIMP
c
      SUBROUTINE ProcessFluxTube(xin,yin,zin,mode,tube_3D)
c      SUBROUTINE ProcessFluxTube(xin,yin,zin,
c     .             mode,tube_3D,
c     .             MAXSURFACE,MAXPOINTS,nsur,npts,hsur,vsur)

      IMPLICIT none

c      INTEGER, INTENT(IN) :: MAXSURFACE,MAXPOINTS
      INTEGER, INTENT(IN) :: mode
      REAL   , INTENT(IN) :: xin,yin,zin
c      INTEGER, INTENT(OUT) :: nsur,npts(MAXSURFACE),hsur(MAXSURFACE)
c      REAL*8 , INTENT(OUT) :: vsur(3,MAXPOINTS,0:MAXSURFACE)

      TYPE    type_tube_3D
         REAL      :: dangle
         INTEGER   :: ncell
         REAL*8    :: area  (1000)
         REAL*8    :: length(1000)
         REAL*8    :: volume(1000)
         INTEGER   :: n
         INTEGER   :: index(  200000)
         REAL*8    :: v    (3,200000)
      ENDTYPE type_tube_3D

      TYPE(type_tube_3D) :: tube_3D(5)

      REAL    FindSeparatrixRadius

      INTEGER i1,i2,itube,ik,ike,nvtx,isegment,nsegments,ring,val_ring
      LOGICAL last_point,polygon_good
      REAL    x1,y1,z1,dradius,rseparatrix,rstart,rend,rvalue(3)   ! *** fix the 1000 and the 10000
      REAL*8  len1,len2,x,z,angle,dangle,vtx(1:3,1000),length(1000),
     .        center(3),vector(3),p1(3),p2(3),p3(3),A,B,C,denominator,u,
     .        dist,mindist,vertex(1:3,5),frac,dot_product,side_length,
     .        cross_product(3),diagonal_42(3),diagonal_53(3),area,
     .        volume,volume_total

      INTEGER n,index(10000)
      REAL    fraction(10000)
      REAL*8  v(3,10000)



      INTEGER ntube_3D
      TYPE(type_tube_3D) :: save_tube_3D(2)


      len1 = 1000.0D0
      len2 = 1000.0D0

      IF (mode.EQ.2) THEN
        save_tube_3D(1) = tube_3D(3)
        save_tube_3D(2) = tube_3D(4)
      ENDIF


      ntube_3D = 0

      dradius = 25.994 * 0.001  ! 0.00054 ! 0.0002
      dangle = 1.0D0

c
c     Looking down the field line:
c
c       4     5
c
c          1
c
c       3     2
c

      rstart = xin !  50.96412 * 0.001
      rend   = yin ! 102.95100 * 0.001
      val_ring = NINT(zin)

      rseparatrix = FindSeparatrixRadius(1)

c      x1 = x1 - dradius
 
      nsegments = 1

      volume_total = 0.0D0

      DO isegment = 1, nsegments

        rvalue(1) = rstart + REAL(isegment-1) / REAL(nsegments) * 
     .                       (rend - rstart)
        rvalue(3) = rstart + REAL(isegment  ) / REAL(nsegments) * 
     .                       (rend - rstart)
        rvalue(2) = 0.5 * (rvalue(1) + rvalue(3))


        WRITE(0,*) 'RVALUES:',rvalue(1:3)
        WRITE(0,*) 'RVALUES:',xin-dradius,xin,xin+dradius
        WRITE(0,*) 'RVALUES:',rvalue(1:3)+rseparatrix
        WRITE(0,*) 'RVALUES:',
     .     rseparatrix+xin-dradius,xin+rseparatrix,
     .     rseparatrix+xin+dradius

        x1 = rseparatrix + rvalue(2)
        y1 = 0.0
        z1 = 0.0

        CALL TraceFieldLine_DIVIMP(x1,y1,z1,2,1,len1,len2,n,v,
     .                             index,fraction,ring,10000)
        IF (ring.NE.val_ring) 
     .    CALL ER('ProcessFluxTube','Invalid ring A',*99)
        ntube_3D = ntube_3D + 1
        tube_3D(ntube_3D)%n = n
        tube_3D(ntube_3D)%index(1:n) = index(1:n)
        DO i1 = 1, n
          tube_3D(ntube_3D)%v(1:3,i1) = v(1:3,i1)
c          IF (index(i1).GT.21) 
c     .      WRITE(0,'(A,2I6,F10.2,2X,3F12.5)')
c     .        'index:',i1,index(i1),fraction(i1),
c     .        SNGL(tube_3D(ntube_3D)%v(1:3,i1))
        ENDDO

c        x1 = x1 + dradius
c        x1 = x1 - dradius
        IF (mode.EQ.1) THEN
          x1 = rseparatrix + rvalue(1)
          CALL TraceFieldLine_DIVIMP(x1,y1,z1,2,1,len1,len2,n,v,
     .                               index,fraction,ring,10000)
          IF (ring.NE.val_ring) 
     .      CALL ER('ProcessFluxTube','Invalid ring B',*99)
          ntube_3D = ntube_3D + 1
          tube_3D(ntube_3D)%n = n
          tube_3D(ntube_3D)%index(1:n) = index(1:n)
          DO i1 = 1, n
          tube_3D(ntube_3D)%v(1:3,i1) = v(1:3,i1)
c          IF (index(i1).GT.21) 
c     .      WRITE(0,'(A,2I6,F10.2,2X,3F12.5)')
c     .        'index:',i1,index(i1),fraction(i1),
c     .        SNGL(tube_3D(ntube_3D)%v(1:3,i1))
          ENDDO
c...      Rotate these points toroidally:
          angle = -0.5D0 * dangle * 3.1415927D0 / 180.0D0
          DO i1 = 1, n
            x = tube_3D(ntube_3D)%v(1,i1)
            z = tube_3D(ntube_3D)%v(3,i1)
            tube_3D(ntube_3D)%v(1,i1) =DCOS(angle) * x - DSIN(angle) * z
            tube_3D(ntube_3D)%v(3,i1) =DSIN(angle) * x + DCOS(angle) * z
          ENDDO
        ELSE
          ntube_3D = ntube_3D + 1
          tube_3D(ntube_3D) = save_tube_3D(1)
        ENDIF

c        x1 = x1 + 2.0 * dradius
        x1 = rseparatrix + rvalue(3)
        CALL TraceFieldLine_DIVIMP(x1,y1,z1,2,1,len1,len2,n,v,
     .                             index,fraction,ring,10000)
        IF (ring.NE.val_ring) 
     .    CALL ER('ProcessFluxTube','Invalid ring C',*99)
        ntube_3D = ntube_3D + 1
        tube_3D(ntube_3D)%n = n
        tube_3D(ntube_3D)%index(1:n) = index(1:n)
        DO i1 = 1, n
          tube_3D(ntube_3D)%v(1:3,i1) = v(1:3,i1)
        ENDDO
c...    Rotate these points toroidally:
        angle = -0.5D0 * dangle * 3.1415927D0 / 180.0D0
        DO i1 = 1, n
          x = tube_3D(ntube_3D)%v(1,i1)
          z = tube_3D(ntube_3D)%v(3,i1)
          tube_3D(ntube_3D)%v(1,i1) = DCOS(angle) * x - DSIN(angle) * z
          tube_3D(ntube_3D)%v(3,i1) = DSIN(angle) * x + DCOS(angle) * z
        ENDDO

        ntube_3D = ntube_3D + 1
        tube_3D(ntube_3D) = tube_3D(ntube_3D-1)      
        angle = dangle * 3.1415927D0 / 180.0D0
        DO i1 = 1, n
          x = tube_3D(ntube_3D)%v(1,i1)
          z = tube_3D(ntube_3D)%v(3,i1)
          tube_3D(ntube_3D)%v(1,i1) = DCOS(angle) * x - DSIN(angle) * z
          tube_3D(ntube_3D)%v(3,i1) = DSIN(angle) * x + DCOS(angle) * z
        ENDDO

        IF (mode.EQ.1) THEN
          ntube_3D = ntube_3D + 1
          tube_3D(ntube_3D) = tube_3D(2)      
          angle = dangle * 3.1415927D0 / 180.0D0
          DO i1 = 1, n
            x = tube_3D(ntube_3D)%v(1,i1)
            z = tube_3D(ntube_3D)%v(3,i1)
            tube_3D(ntube_3D)%v(1,i1) =DCOS(angle) * x - DSIN(angle) * z
            tube_3D(ntube_3D)%v(3,i1) =DSIN(angle) * x + DCOS(angle) * z
          ENDDO
        ELSE
          ntube_3D = ntube_3D + 1
          tube_3D(ntube_3D) = save_tube_3D(2)
        ENDIF


c        IF (mode.EQ.1) THEN
c          DO i2 = 1, ntube_3D
c            IF (mode.EQ.2.and.(i1.EQ.2.OR.i1.EQ.5)) CYCLE
c            DO i1 = 1, tube_3D(i2)%n-1
c              nsur = nsur + 1
c              hsur(nsur) = -1
c              npts(nsur) =  2
c              vsur(1:3,1,nsur) = tube_3D(i2)%v(1:3,i1  )
c              vsur(1:3,2,nsur) = tube_3D(i2)%v(1:3,i1+1)
c            ENDDO
c          ENDDO
c        ENDIF

c       CYCLE

c      RETURN




       itube = 1
       ike = 0
       DO i1 = 1, tube_3D(itube)%n 
         ike = MAX(ike,tube_3D(itube)%index(i1))
       ENDDO
       WRITE(0,*) 'IKE:',ike

       tube_3D(1)%ncell = ike
       tube_3D(1)%dangle = dangle
       tube_3D(1)%area  (1:ike) = 0.0
       tube_3D(1)%length(1:ike) = 0.0
       tube_3D(1)%volume(1:ike) = 0.0

c...  START OF LOOP:
c        DO ik = 20, 36
        DO ik = 1, ike

c...      Collect line segments associated with this cell:
          last_point = .FALSE.
          nvtx = 0

        
c...      Select the line segments that are associated with the
c         cell IK:
          itube = 1
          DO i1 = 1, tube_3D(itube)%n 
            IF (tube_3D(itube)%index(i1).EQ.ik) THEN
              last_point = .TRUE.
              nvtx = nvtx + 1
              vtx(1:3,nvtx) =  tube_3D(itube)%v(1:3,i1)
            ELSEIF (last_point.AND.ik.LT.ike) THEN
              nvtx = nvtx + 1
              vtx(1:3,nvtx) = tube_3D(itube)%v(1:3,i1)          
              EXIT
            ENDIF
          ENDDO
        
          DO i1 = 1, nvtx  
c            WRITE(0,'(A,I6,3F12.5)') 'VTX:',i1,SNGL(vtx(1:3,i1))
          ENDDO

c...      Calulate the total length:
          length = 0.0D0
          DO i1 = 1, nvtx-1
            length(i1+1) = length(i1) +
     .               DSQRT((vtx(1,i1+1)-vtx(1,i1))**2 +
     .                     (vtx(2,i1+1)-vtx(2,i1))**2 +
     .                     (vtx(3,i1+1)-vtx(3,i1))**2)
          ENDDO
c...      Find half way point:
          DO i2 = 1, nvtx-1
            IF (length(i2  ).LE.0.5D0*length(nvtx).AND.
     .          length(i2+1).GT.0.5D0*length(nvtx)) THEN
              frac = (0.5D0*length(nvtx) - length(i2)) / 
     .               (      length(i2+1) - length(i2))
              center(1:3) = (1.0D0 - frac) * vtx(1:3,i2  ) +
     .                               frac  * vtx(1:3,i2+1)
              vector(1:3) = vtx(1:3,i2+1) - vtx(1:3,i2)
              EXIT
            ENDIF
          ENDDO
c          DO i1 = 1, nvtx
c            WRITE(0,*) 'LENGTH:',i1,length(i1),i2,frac
c          ENDDO
          WRITE(0,*) 'CENTER:',center(1:3)
          WRITE(0,*) 'VECTOR:',vector(1:3)

          vector(1:3) = vector(1:3) /
     .                 DSQRT(vector(1)**2 + vector(2)**2 + vector(3)**2)
        
c...      Find intersections between the flux-tube boundary field line tracings
c         and the plane defined by the center point of the cell:
          itube = 2
          A = vector(1)
          B = vector(2)
          C = vector(3)
          vertex = 0.0D0
          vertex(1:3,1) = center(1:3)
          DO itube = 2, 5
            mindist = 1.0D+20
            DO i1 = 1, tube_3D(itube)%n-1
              p1(1:3) = tube_3D(itube)%v(1:3,i1  ) - center(1:3)
              p2(1:3) = tube_3D(itube)%v(1:3,i1+1) - center(1:3)
              denominator = A * (p1(1) - p2(1)) + B * (p1(2) - p2(2)) + 
     .                      C * (p1(3) - p2(3))
              IF (DABS(denominator).LT.1.0D-10) THEN
                u = -1.0D0
              ELSE
                u = (A * p1(1) + B * p1(2) + C * p1(3)) / denominator
              ENDIF
              IF (u-1.0D-10.GT.0.0D0.AND.u+1.0D-10.LT.1.0D0) THEN
                p3(1:3) = p1(1:3) + u * (p2(1:3) - p1(1:3))
                dist = DSQRT(p3(1)**2 + p3(2)**2 + p3(3)**2)
                IF (dist.LT.mindist) THEN
                  vertex(1:3,itube) = p3(1:3) + center(1:3)
                  mindist = dist
c                  WRITE(0,*) 'U:',i1,u,dist
c                  WRITE(0,*) 'P3:',p3(1:3) + center(1:3)
c                  WRITE(0,*) 'BINGO!',itube
                ENDIF
              ENDIF
            ENDDO
          ENDDO

c...      Check that the plane made up of the identified points is indeed
c         perpendicular to the central vector of the cell:
          DO i1 = 1, 4
            i2 = i1 + 1
            IF (i2.EQ.5) i2 = 1        
            dot_product = 
     .        vector(1) * (vertex(1,i1)-vertex(1,i2)) + 
     .        vector(2) * (vertex(2,i1)-vertex(2,i2)) + 
     .        vector(3) * (vertex(3,i1)-vertex(3,i2))
            IF (dot_product.GT.1.0D-10) THEN
              WRITE(0,*) 'CROSS_PRODUCT PAIN'
              STOP
            ENDIF
          ENDDO

c...      Check the length of each side:
          polygon_good = .TRUE.
          DO i1 = 2, 5
            i2 = i1 + 1
            IF (i2.EQ.6) i2 = 1               
            side_length = DSQRT((vertex(1,i1)-vertex(1,i2))**2 + 
     .                          (vertex(2,i1)-vertex(2,i2))**2 + 
     .                          (vertex(3,i1)-vertex(3,i2))**2)
            WRITE(0,*) 'SIDE LENGTH:',i1-1,side_length
            IF (side_length.GT.0.5) polygon_good = .FALSE.
          ENDDO

          WRITE(0,*) 'POLYGON_GOOD:',polygon_good

          IF (.FALSE.) THEN
            length(nvtx) = 1.0

            vertex(1,2) = 1.0
            vertex(2,2) = 0.0
            vertex(3,2) = 0.0

            vertex(1,3) = 1.0
            vertex(2,3) = 1.0
            vertex(3,3) = 0.0

            vertex(1,4) = 2.0
            vertex(2,4) = 1.0
            vertex(3,4) = 0.0

            vertex(1,5) = 2.0
            vertex(2,5) = 0.0
            vertex(3,5) = 0.0
 
            vector(1) = 0.0
            vector(2) = 0.0
            vector(3) = 1.0
          ENDIF

          IF (polygon_good) THEN
! AREA of a 3D quadrelateral, n is unit normal vector of the plane...
      ! 2 A = n DOT ( v2 - v0 cross v3 - v1 ) 
          
c...       Calculate the area of a quadrelateral:
           diagonal_42(1:3) = vertex(1:3,4) - vertex(1:3,2)
           diagonal_53(1:3) = vertex(1:3,5) - vertex(1:3,3)           
        
           cross_product(1) = diagonal_42(2) * diagonal_53(3) -   ! I had a function for this somewhere...
     .                        diagonal_42(3) * diagonal_53(2)
           cross_product(2) = diagonal_42(3) * diagonal_53(1) -
     .                        diagonal_42(1) * diagonal_53(3)
           cross_product(3) = diagonal_42(1) * diagonal_53(2) -
     .                        diagonal_42(2) * diagonal_53(1)
        
           area = DABS(0.5D0 * (vector(1) * cross_product(1) + 
     .                          vector(2) * cross_product(2) + 
     .                          vector(3) * cross_product(3)))

          ELSE

            area = 0.0D0

          ENDIF

          volume = area * length(nvtx)

          volume_total = volume_total + volume

          WRITE(0,*) 'AREA,VOLUME:',ik,area,
     .               volume*(360.0D0/dangle)  ! *** LEFT OFF ***

          tube_3D(1)%area  (ik) = area
          tube_3D(1)%length(ik) = length(nvtx)
          tube_3D(1)%volume(ik) = volume

c          DO i1 = 2, 5
c            i2 = i1 + 1
c            IF (i2.EQ.6) i2 = 2
c            nsur = nsur + 1
c            hsur(nsur) = -2
c            npts(nsur) =  2
c            vsur(1:3,1,nsur) = vertex(1:3,i1)
c            vsur(1:3,2,nsur) = vertex(1:3,i2)
c          ENDDO
    
        ENDDO

      ENDDO

c      WRITE(0,*) '*** TOTAL VOLUME ***',volume_total*(360.0D0/dangle)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: SelectGridRegoin_DIVIMP
c
      SUBROUTINE SelectGridRegion_DIVIMP(rhoval,nrings,rings,MAX_IR)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER, INTENT(IN)  :: MAX_IR
      REAL   , INTENT(IN)  :: rhoval
      INTEGER, INTENT(OUT) :: nrings,rings(MAX_IR)

      INTEGER ik,ir,ir1,i1,i2
      LOGICAL addtolist

c...  Find which ring the particular value of RHO corresponds to:
      DO ir = 2, irwall
        IF (rho(ir,CELL1).EQ.0.0) CYCLE
c...    Take the first case since RHO in some secondary PFR could catch things up:
        IF (rhoval.GE.rho(ir,IN14).AND.rhoval.LT.rho(ir,OUT23)) EXIT  
      ENDDO
      IF (ir.EQ.irwall+1) 
     .  CALL ER('SelectGridRegion_DIVIMP','Ring not identified',*99)

      nrings = 1
      rings(nrings) = ir

c...  Assemble a list of immediate neighbours:
      DO ik = 1, nks(ir)
        DO i1 = 1, 2
          IF (i1.EQ.1) ir1 = irins (ik,ir)
          IF (i1.EQ.2) ir1 = irouts(ik,ir)
          addtolist = .TRUE.
          DO i2 = 1, nrings
            IF (rings(i2).EQ.ir1) addtolist = .FALSE.
          ENDDO
          IF (addtolist) THEN
            nrings = nrings + 1
            IF (nrings.GT.MAX_IR) 
     .        CALL ER('SelectGridRegion_DIVIMP','IR bound exceeded',*99)
            rings(nrings) = ir1
          ENDIF
        ENDDO
      ENDDO

      RETURN
99    WRITE(0,*) '  RHOVAL= ',rhoval
      WRITE(0,*) '  RHO   = ',rho(2:irsep,CELL1)
      CALL OutputData(85,'FRUSTRATION...')
      STOP
      END
c
c ======================================================================
c
      REAL FUNCTION FindSeparatrixRadius(mode)   
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER, INTENT(IN) :: mode

      INTEGER FindMidplaneCell

      INTEGER ikm,ik,ir,id
      REAL    rmid

c...  Check if the grid is connected double-null -- lame, should
c     really pass the value of CONNECTED from DIVIMP:
      IF (rs(ikto,irsep).LT.rxp.AND.rs(ikti,irsep).LT.rxp) 
     .  connected = .TRUE.

      IF (connected) THEN
        ir = irsep2
      ELSE
        ir = irsep
      ENDIF

      ikm = FindMidplaneCell(ir)
      id = korpg(ikm,ir)
      rmid = MAX(rvertp(1,id),rvertp(4,id))

      IF (.NOT.connected.AND.rmid.LT.rxp) THEN
c...    Really want the outer midplane radius, but this result
c       suggests the grid is connected, so try again:
        ikm = FindMidplaneCell(irsep2)
        id = korpg(ikm,irsep2)
        rmid = MAX(rvertp(1,id),rvertp(4,id))
        WRITE(0,*) 'WHOA! Looks like a connected grid but CONNECTED '//
     .             'not set'
        STOP 'HALTING CODE'
      ENDIF

      FindSeparatrixRadius = rmid

      RETURN
99    STOP
      END
c
c ======================================================================
c
      INTEGER FUNCTION FindMidplaneCell(ir)   
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'

      INTEGER, INTENT(IN) :: ir

      INTEGER ikm,ik,id
      REAL    rmid,z(2)

c...  Find midplane cell:   ! *** NEED TO DO THIS SO THAT PHI IS REFERENCED PROPERLY TO PHI=0 AT THE OUTER MIDPLANE ***
      ikm = -1              !     NEED A BETTER WAY...
      rmid = -100.0
      DO ik = 1, nks(ir)
        id = korpg(ik,ir)
        z(1) = 0.5 * (zvertp(1,id) + zvertp(2,id))
        z(2) = 0.5 * (zvertp(3,id) + zvertp(4,id))
        IF (((z(1).LT.z0.AND.z(2).GE.z0).OR.
     .       (z(2).LT.z0.AND.z(1).GE.z0)).AND.
c        IF (((z(1).LT.0.0.AND.z(2).GE.0.0).OR.
c     .       (z(2).LT.0.0.AND.z(1).GE.0.0)).AND.
     .      rs(ik,ir).GT.rmid) THEN
          ikm  = ik
          rmid = rs(ik,ir)
        ENDIF
      ENDDO

      FindMidplaneCell = ikm

      RETURN
99    STOP
      END
c
c
c
c ======================================================================
c SHOULD NOT BE HERE!  TEMPORARY...
c ======================================================================
c
c *** AM I USING BRRATIO CORRECTLY? DO I NEED TO MAP THE ANGLE SOMEHOW, 
c SO THAT I'M GETTING THE PITCH ANGLE RIGHT IN THE FIELD LINE
c COORDINATE SYSTEM? IS THIS ALL BOGUS? ***
c
c subroutine: TraceFieldLine_DIVIMP
c
      SUBROUTINE TraceFieldLine_DIVIMP(xin,yin,zin,mode,chop,
     .                                 length1,length2,
     .                                 n,v,index,fraction,ring,MAXN)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'colours'
      INCLUDE 'slout'
      INCLUDE 'slcom'

      INTEGER FindMidplaneCell
      REAL    FindSeparatrixRadius

      INTEGER, INTENT(IN) :: mode,MAXN,chop
      REAL   , INTENT(IN) :: xin,yin,zin
      REAL*8              :: length1,length2
c      INTEGER, INTENT(IN) :: region,mode,MAXN
c      REAL   , INTENT(IN) :: rad_coord,in_phi
      INTEGER, INTENT(OUT) :: n,index(MAXN),ring
      REAL   , INTENT(OUT) :: fraction(MAXN)
      REAL*8 , INTENT(OUT) :: v(3,MAXN)


      INTEGER fp,ik,ir,ikm,id,iobj,isur,ipts,i1,ik1,ir1,ike
      LOGICAL finished
      REAL    rsep
      REAL*8  r(2),z(2),deltar,deltaz,deltac,deltap,phi,dphi,rval,pval,
     .        frac1,frac2,angle1,angle2,brat,rfrac,r12frac,z12frac,
     .        rhoin,phiin,frac,rpos,rvalmin,zvalmin,zvalmax,
     .        len1,len2,lenmax1,lenmax2,
     .        dpol,dtor,alpha

      REAL       TOL
      PARAMETER (TOL=1.0D-06)


      fp = 0

      dphi = 0.5 !  1.0 ! 5.0 ! 10.0  ! Make this an adjustable parameter...

      
      IF (yin.NE.0.0) 
     .  CALL ER('TraceFieldLine_DIVIMP','Sorry, midplane only',*99)

c...  Convert from Cartesean coordinates to RHO and PHI:
      rsep = FindSeparatrixRadius(1)
      rhoin = DBLE(SQRT(xin**2 + zin**2) - rsep)
      IF (ABS(xin).LT.1.0E-06) THEN
        IF (zin.GT.0.0) THEN
          phiin = 90.0
        ELSE
          phiin = 270.0
        ENDIF
      ELSE
        phiin = DBLE(ATAN(ABS(zin / xin)) * 180.0 / PI)
        IF (xin.LT.0.0.and.zin.GT.0.0) phiin = 180.0 - phiin
        IF (xin.LT.0.0.and.zin.LT.0.0) phiin = 180.0 + phiin
        IF (xin.GT.0.0.and.zin.LT.0.0) phiin = 360.0 - phiin
      ENDIF
c      rhoin = xin  ! *** TEMP ***
c      phiin = zin  ! *** TEMP ***
c      WRITE(0,*) '==INPUT:',rhoin,phiin
c      WRITE(0,*) '       :',xin,zin,rsep
c      STOP 'sdfsdf'

      n = 0

      len1 = 0.0D0
      len2 = 0.0D0

      SELECTCASE (chop)
        CASE(1)    ! No restrictions
          rvalmin =  0.0D0
          zvalmin = -1.0D+20
          zvalmax =  1.0D+20
          lenmax1 =  1.0D+20
          lenmax2 =  1.0D+20
        CASE(2:3,6)  ! No field lines above and below or inside the low field side, main plasma separatrix 
          rvalmin =  0.0D0
          zvalmin = -1.0D+20
          zvalmax =  1.0D+20
          ir = irsep-1
          IF (idring(ir).NE.BOUNDARY) THEN
            rvalmin = DBLE(rxp)
            zvalmin =  1.0D+20
            zvalmax = -1.0D+20
            DO ik = 1, nks(ir)-1
              DO i1 = 1, 4
                id = korpg(ik,ir)
                zvalmin = MIN(zvalmin,DBLE(zvertp(i1,id)))
                zvalmax = MAX(zvalmax,DBLE(zvertp(i1,id)))
              ENDDO
            ENDDO
          ENDIF
          IF (chop.EQ.6) THEN
            lenmax1 = length1
            lenmax2 = length2
          ELSE
            lenmax1 = 1.0D+20
            lenmax2 = 1.0D+20
          ENDIF
c          WRITE(0,*) 'ZVALs:',zvalmin,zvalmax
        CASE(4:5)  ! Restricted length...
          rvalmin =  0.0D0
          zvalmin = -1.0D+20
          zvalmax =  1.0D+20
          lenmax1 =  length1
          lenmax2 =  length2
        CASE DEFAULT
          CALL ER('TraceFieldLine_DIVIMP','Unrecognised CHOP',*99)
      ENDSELECT
c      WRITE(0,*) 'CHOP:',chop,rvalmin,zvalmin,zvalmax

c...  Find where we are on the outer midplane:
      DO ir = 2, irwall

c...    For now, ignore the inner SOL for a CDN grid (CONNECTED
c       currently assigned/hacked in FindSeparatrixRadius, which
c       has to be called before this check is made, nasty):
        IF (connected.AND.ir.GE.irsep.AND.ir.LT.irsep2) CYCLE

        IF (rho(ir,CELL1).EQ.0.0) CYCLE

        IF (rhoin.GE.DBLE(rho(ir,SIDE14)).AND.
     .      rhoin.LT.DBLE(rho(ir,SIDE23))) THEN
c...      Find midplane cell:   ! *** NEED TO DO THIS SO THAT PHI IS REFERENCED PROPERLY TO PHI=0 AT THE OUTER MIDPLANE ***
          ikm = -1              !     NEED A BETTER WAY...
          ikm = FindMidplaneCell(ir)
          IF (ikm.EQ.-1) CALL ER('TraceFieldLine_DIVIMP','No '//
     .                           'midplane cell found',*99)
          
c...      Decide how to assign the magnetic field line pitch angle information:
          SELECTCASE (mode)
            CASE (1)
c...          Just the pitch angle at the center of the cell all the time:
              ir1 = 0
              frac = DBLE(rhoin          - rho(ir,SIDE14)) / 
     .               DBLE(rho(ir,SIDE23) - rho(ir,SIDE14))
c              frac = 0.5D0
              r12frac = frac
              z12frac = frac ! 0.5D0
              rfrac = 0.0D0
            CASE (2)
c...          Interpolate the field line pitch angle, gives a continuous :
              frac = DBLE(rhoin          - rho(ir,SIDE14)) / 
     .               DBLE(rho(ir,SIDE23) - rho(ir,SIDE14))
              r12frac = frac
              z12frac = frac
              IF (rhoin.LT.rho(ir,CELL1)) THEN
                ir1 = irins(ikm,ir)

                WRITE(0,*) 'IRs:',ir,ir1
                rfrac =     -(rhoin          - DBLE(rho(ir,CELL1))) / 
     .                   DBLE(rho(ir1,CELL1) -      rho(ir,CELL1) )
              ELSE
                ir1 = irouts(ikm,ir)
                rfrac =      (rhoin          - DBLE(rho(ir,CELL1))) / 
     .                   DBLE(rho(ir1,CELL1) -      rho(ir,CELL1) )
c                WRITE(fp,*) rhoin,rho(ir,CELL1),rho(ir1,CELL1)
              ENDIF
            CASE DEFAULT
              CALL ER('TraceFieldLine_DIVIMP','Unrecognised MODE',*99)
          ENDSELECT

          WRITE(fp,'(A,2I6,3F10.4)') 
     .      ' ==MIDPLANE:',ikm,ir,REAL(frac),REAL(rfrac),REAL(rhoin)
          EXIT
        ENDIF
      ENDDO

      ring = ir

c...  Work from midplane to low IK target:
      phi = phiin
      DO ik = ikm, 1, -1  ! 1, -1
        id = korpg(ik,ir)
        r(1) =          r12frac  * DBLE(rvertp(3,id)) + 
     .         (1.0D0 - r12frac) * DBLE(rvertp(4,id))
        z(1) =          z12frac  * DBLE(zvertp(3,id)) + 
     .         (1.0D0 - z12frac) * DBLE(zvertp(4,id))
        r(2) =          r12frac  * DBLE(rvertp(2,id)) +
     .         (1.0D0 - r12frac) * DBLE(rvertp(1,id))
        z(2) =          z12frac  * DBLE(zvertp(2,id)) + 
     .         (1.0D0 - z12frac) * DBLE(zvertp(1,id))
        rpos = 0.5D0 * (r(1) + r(2))
        IF (rfrac.LE.0.0D0) THEN
          ik1 = ikins(ik,ir)
          ir1 = irins(ik,ir)
c          WRITE(0,*) 'IK,IR,IK1,IR1=',ik,ir,ik1,ir1
          brat = (1.0D0 + rfrac) * DBLE(bratio(ik,ir)  ) - 
     .                    rfrac  * DBLE(bratio(ik1,ir1))
        ELSE
          ik1 = ikouts(ik,ir)
          ir1 = irouts(ik,ir)
          brat = (1.0D0 - rfrac) * DBLE(bratio(ik,ir)  ) + 
     .                    rfrac *  DBLE(bratio(ik1,ir1))

        ENDIF
c        WRITE(0,'(A,5F10.5)') 
c     .    'BRAT:',REAL(brat),REAL(r(1:2)),REAL(z(1:2))
        deltar = r(2) - r(1)
        deltaz = z(2) - z(1)
        dpol = DSQRT(deltar**2 + deltaz**2)
        alpha = DASIN(brat)
        dtor = dpol / DTAN(alpha)
        deltac = dtor
c        deltac = ABS(deltaz) / brat
        deltap = -1.0D0 * deltac / rpos * 180.0D0 / DBLE(PI)
        angle1 = 0.0D0
        finished = .FALSE.
        DO WHILE (angle1.GT.deltap)
          angle2 = MAX(angle1-dphi,deltap) 
          frac1 = angle1 / deltap
          frac2 = angle2 / deltap
          IF (angle1.EQ.0.0D0.AND.ik.EQ.ikm) THEN
            n = n + 1
            IF (n.GT.MAXN) CALL ER('TraceFieldLine_DIVIMP','N bust',*99)
            v(1,n) = r(1)
            v(2,n) = z(1)
            v(3,n) = phi
            index(n) = ik + 1
            fraction(n) = frac2
c...        Convert from r,z,phi to x,y,z (y okay already):
            rval = v(1,n)
            pval = v(3,n) * DBLE(PI) / 180.0D0
            v(1,n) = rval * DCOS(pval)
            v(3,n) = rval * DSIN(pval)
c            WRITE(0,*) '==START OF FILAMENT:',SNGL(v(3,n)),SNGL(phi)
          ENDIF
c...      Don't follow the field line at all if length restriction set to zero:
          IF ((chop.EQ.4.OR.chop.EQ.5.OR.chop.EQ.6).AND.
     .        lenmax1.EQ.0.0D0) THEN
            finished = .TRUE.
            EXIT
          ENDIF
          n = n + 1
          IF (n.GT.MAXN) CALL ER('TraceFieldLine_DIVIMP','N bust',*99)
          v(1,n) = r(1) + frac2 * deltar
          v(2,n) = z(1) + frac2 * deltaz
          v(3,n) = phi + angle2
          index(n) = ik
          fraction(n) = frac2
c...      Convert from r,z,phi to x,y,z (y okay already):
          rval = v(1,n)
          pval = v(3,n) * DBLE(PI) / 180.0D0
          v(1,n) = rval * DCOS(pval)
          v(3,n) = rval * DSIN(pval)
          len1 = len1 + DSQRT((v(1,n)-v(1,n-1))**2 +
     .                        (v(2,n)-v(2,n-1))**2 +
     .                        (v(3,n)-v(3,n-1))**2)

c          WRITE(0,*) 'ZVAL-:',v(2,n),zvalmax,zxp
c          WRITE(0,*) 'ZVAL1-:',ik,frac2

          IF ((rval  .LE.rvalmin).OR.
     .        (v(2,n).LE.zvalmin.OR.v(2,n).GE.zvalmax).OR.
c     .        (zxp.LT.0.0.AND.v(2,n).LE.zvalmin).OR.
c     .        (zxp.GT.0.0.AND.v(2,n).GE.zvalmax).OR.
     .        (len1  .GE.lenmax1)) THEN
            finished = .TRUE.
            EXIT
          ENDIF
          angle1 = angle2
        ENDDO 
        IF (finished) EXIT
        phi = phi + deltap
      ENDDO

c...  Swap order of these points, so that they start at the low IK target
c     and proceed to the midplane:
      DO i1 = 1, n/2
        v(1:3,n+1   ) = v(1:3,i1    )
        v(1:3,i1    ) = v(1:3,n-i1+1)
        v(1:3,n-i1+1) = v(1:3,n+1   )
        index(n+1   ) = index(i1    )
        index(i1    ) = index(n-i1+1)
        index(n-i1+1) = index(n+1   )
        fraction(n+1   ) = fraction(i1    )
        fraction(i1    ) = fraction(n-i1+1)
        fraction(n-i1+1) = fraction(n+1   )
      ENDDO

c...  Work from midplane to high IK target:
      phi = phiin
      ike = nks(ir)
      IF (ir.LT.irsep) ike = ike - 1
      DO ik = ikm+1, ike
c...    Don't follow the field line at all if length restriction set to zero:
        IF ((chop.EQ.4.OR.chop.EQ.5.OR.chop.EQ.6).AND.
     .      lenmax2.EQ.0.0D0) EXIT
        id = korpg(ik,ir)
c        frac = 0.5D0 * (1.0D0 + rfrac)
        r(1) =          r12frac  * DBLE(rvertp(2,id)) + 
     .         (1.0D0 - r12frac) * DBLE(rvertp(1,id))
        z(1) =          z12frac  * DBLE(zvertp(2,id)) + 
     .         (1.0D0 - z12frac) * DBLE(zvertp(1,id))
        r(2) =          r12frac  * DBLE(rvertp(3,id)) +
     .         (1.0D0 - r12frac) * DBLE(rvertp(4,id))
        z(2) =          z12frac  * DBLE(zvertp(3,id)) + 
     .         (1.0D0 - z12frac) * DBLE(zvertp(4,id))
        rpos = 0.5D0 * (r(1) + r(2))
        IF (rfrac.LE.0.0D0) THEN
          ik1 = ikins(ik,ir)
          ir1 = irins(ik,ir)
          brat = (1.0D0 + rfrac) * DBLE(bratio(ik,ir  )) -
     .                    rfrac  * DBLE(bratio(ik1,ir1))
        ELSE
          ik1 = ikouts(ik,ir)
          ir1 = irouts(ik,ir)
          brat = (1.0D0 - rfrac) * DBLE(bratio(ik,ir  )) + 
     .                    rfrac  * DBLE(bratio(ik1,ir1))
        ENDIF
c        brat = DBLE(bratio(ik,ir))
        deltar = r(2) - r(1)
        deltaz = z(2) - z(1)
        dpol = DSQRT(deltar**2 + deltaz**2)
        alpha = DASIN(brat)
        dtor = dpol / DTAN(alpha)
        deltac = dtor
c        deltac = ABS(deltaz) / brat
        deltap = deltac / rpos * 180.0D0 / DBLE(PI)
        angle1 = 0.0D0
        finished = .FALSE.
        DO WHILE (angle1.LT.deltap)
          angle2 = MIN(angle1+dphi,deltap) 
          frac1 = angle1 / deltap
          frac2 = angle2 / deltap
          n = n + 1
          IF (n.GT.MAXN) CALL ER('TraceFieldLine_DIVIMP','N bust',*99)
          v(1,n) = r(1) + frac2 * deltar
          v(2,n) = z(1) + frac2 * deltaz
          v(3,n) = phi + angle2
          index(n) = ik
          fraction(n) = frac2
c...      Convert from r,z,phi to x,y,z (y okay already):
          rval = v(1,n)
          pval = v(3,n) * DBLE(PI) / 180.0D0
          v(1,n) = rval * DCOS(pval)
          v(3,n) = rval * DSIN(pval)
          len2 = len2 + DSQRT((v(1,n)-v(1,n-1))**2 +
     .                        (v(2,n)-v(2,n-1))**2 +
     .                        (v(3,n)-v(3,n-1))**2)
          IF ((rval  .LE.rvalmin).OR.
     .        (v(2,n).LE.zvalmin.OR.v(2,n).GE.zvalmax).OR.
c     .        (zxp.LT.0.0.AND.v(2,n).GE.zvalmax).OR.
c     .        (zxp.GT.0.0.AND.v(2,n).LE.zvalmin).OR.
     .        (len2  .GE.lenmax2 )) THEN
            finished = .TRUE.
            EXIT
          ENDIF

c          WRITE(0,*) 'ZVAL2-:',ik,frac2,ike

          angle1 = angle2
        ENDDO 
        IF (ik.LT.ike) index(n) = ik + 1
        IF (finished) EXIT
        phi = phi + deltap
      ENDDO

c      DO i1 = 1, n
c        WRITE(0,*) '-->',i1,ike,index(i1)
c      ENDDO

c...  Convert from r,z,phi to x,y,z (y okay already):
c      DO i1 = 1, n
c        rval = v(1,i1)
c        pval = v(3,i1) * DBLE(PI) / 180.0D0
c        v(1,i1) = rval * DSIN(pval)
c        v(3,i1) = rval * DCOS(pval)
cc        WRITE(fp,*) ' V:',v(1:3,i1)
c      ENDDO

c      WRITE(0,*) 'n:',n

c      WRITE(0,*) 'LENGTH:',len1,len2

      SELECTCASE (chop)
        CASE(1)
          length1 = len1
          length2 = len2
        CASE(2)
        CASE(3)
          length1 = len1
          length2 = len2
        CASE(4)
        CASE(5)
          length1 = len1
          length2 = len2
        CASE(6)
        CASE DEFAULT
          CALL ER('TraceFieldLine_DIVIMP','Unrecognised CHOP',*99)
      ENDSELECT

      RETURN
 99   WRITE(0,*) ' CHOP= ',chop
      STOP
      END
