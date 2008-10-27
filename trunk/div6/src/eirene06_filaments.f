c
c ======================================================================
c
c
c *** Need to look over zones maybe since each fully split tetrahedron
c     goes into 24 tetrahedrons, which may blow the gasket *** 
c
c
      SUBROUTINE ResolveFilament
      USE mod_geometry
      USE mod_filament_params
      USE mod_filament
      USE mod_eirene06_locals
      USE mod_eirene06  ! *** TEMP ***
      IMPLICIT none

c      REAL FindSeparatrixRadius

      INTEGER iobj,iobj1,tmp_nobj,ifilament,ivtx,icell,iside,mode,iloop

      INTEGER n,nir,nmid
      REAL    x,y,z,rho,phi
      REAL*8  len1,len2,v(3,10000),scale
      

      obj(1:nobj)%segment(1) = 0

c      RETURN

      SELECTCASE (3)
        CASE(-1)
        CASE(0)
          DO iobj = 1, nobj
            IF (obj(iobj)%index(IND_IK).NE.15.OR.
     .          (obj(iobj)%index(IND_IR).NE.5.AND.
     .           obj(iobj)%index(IND_IR).NE.6).OR.
     .          obj(iobj)%index(IND_IS).NE.3) CYCLE
            obj(iobj)%segment(1) = 1
          ENDDO
          tmp_nobj = nobj
          DO iobj = 1, tmp_nobj
            IF (obj(iobj)%segment(1).GT.0) CALL DivideTetrahedron(iobj)
          ENDDO
          CALL CleanObjectArray

        CASE(1)
          DO iobj = 1, nobj
            IF (obj(iobj)%index(IND_IS).EQ.3) obj(iobj)%segment(1) = 1
            IF (obj(iobj)%index(IND_IS).EQ.4) obj(iobj)%segment(1) = 1
          ENDDO
          tmp_nobj = nobj
          DO iobj = 1, tmp_nobj
            IF (obj(iobj)%segment(1).GT.0) CALL DivideTetrahedron(iobj)
          ENDDO
          CALL CleanObjectArray

        CASE (2)
          CALL TraceFieldLine_DIVIMP(0.0,0.0,0.0,1,1,1.0E+20,1.0E+20,
     .                               n,v,10000)
          CALL SelectTetrahedrons(n,v(1:3,1:n))
          tmp_nobj = nobj
          DO iobj = 1, tmp_nobj
            IF (obj(iobj)%segment(1).GT.0) CALL DivideTetrahedron(iobj)
          ENDDO
          CALL CleanObjectArray

        CASE (3)
c          CALL DefineFilaments 
c          CALL SetupFilaments(0.0D0)

c          rsep = FindSeparatrixRadius(1)

          DO iloop = 1, 2
            obj(1:nobj)%segment(1) = 0

            mode = 1              ! Just compare perpendicular distance between 
            IF (iloop.EQ.1) THEN  ! line segment and tetrahedron centroid 
              scale = 0.05D0      ! For MODE=1, distance threshold for selection
            ELSE
              scale = 0.02D0
            ENDIF
            nir = 1
c             nir = filament(ifilament)%ir_space(0,icell)

            DO ifilament = 1, nfilament
              DO icell = 1, filament(ifilament)%ncell
c...            Once cells are setup properly, should really have coarser cells stored in FILAMENT% and then temporarily 
c               subdivide them when looking for tetrahedron intersections, rather than always storing so many vertices in 
c               filament cross-section...
                ivtx = icell  ! *** HACK ***
                x = SNGL(filament(ifilament)%vtx(1,ivtx))
                y = 0.0
                z = SNGL(filament(ifilament)%vtx(3,ivtx))
                len1 = filament(ifilament)%lcell(1,ivtx) 
                len2 = filament(ifilament)%lcell(2,ivtx) 
                WRITE(0,*) '==FILAMENT:',ifilament,ivtx,x,z
                CALL TraceFieldLine_DIVIMP(x,y,z,2,4,len1,len2,
     .                                     n,v,10000)  ! *** HACK *** (the 10000)
c                CALL TraceFieldLine_DIVIMP(rho,0.0,phi,2,2,1.0,1.0,
c     .                                     n,v,10000)  ! *** HACK *** (the 10000)
c                CALL TraceFieldLine_DIVIMP(x,y,z,2,2,1.0E+20,1.0E+20,
c     .                                     n,v,10000)  ! *** HACK *** (the 10000)
                CALL SelectTetrahedrons(n,v(1:3,1:n),mode,scale,
     .                 nir,filament(ifilament)%ir_space(1:nir,icell))
              ENDDO  ! TUBE
            ENDDO  ! FILAMENT

            IF (iloop.EQ.1) THEN
c...          Resolve the nighbouring tetrahedrons as well to make sure that the
c             tetrahedrons along the filament are fully resolved (or as resolved as
c             possible anyway):
              DO iobj = 1, nobj                       
                IF (obj(iobj)%segment(1).NE.1) CYCLE  
                DO iside = 1, obj(iobj)%nside         
                 iobj1 = obj(iobj)%omap(iside)
                 IF (iobj1.LE.0) CYCLE
                 IF (obj(iobj1)%segment(1).EQ.0) 
     .             obj(iobj1)%segment(1) = 2  ! *** HACK ***
                ENDDO
              ENDDO           
c...          Increase spatial resolution along the filament:
              tmp_nobj = nobj
              DO iobj = 1, tmp_nobj
                IF (obj(iobj)%segment(1).GT.0) 
     .            CALL DivideTetrahedron(iobj)
              ENDDO
            ENDIF

            CALL CleanObjectArray
          ENDDO  ! LOOP

      ENDSELECT



c     Doesn't work here for some reason, so hack is in ProcessFluidGrid_06...     
c      cell(1)%plasma(1) = 1.0                   ! Te (eV)
c      cell(1)%plasma(2) = 1.0                   ! Ti (eV)
c      cell(1)%plasma(3) = 5.0E+21               ! ni (eV) (ne=ni assumed at present)
      DO iobj = 1, nobj
        IF (obj(iobj)%segment(1).EQ.1) obj(iobj)%index(IND_PLASMA) = 1
      ENDDO



      RETURN
99    STOP
      END
c
c ======================================================================
c ======================================================================
c
c subroutine: SelectTetrahedrons
c
c
      SUBROUTINE SelectTetrahedrons(n,v,mode,scale,n_ir,ir)
      USE mod_geometry
      USE mod_filament
      USE mod_eirene06_locals
      IMPLICIT none

      REAL*8 CalcPerp

      INTEGER, INTENT(IN) :: n,mode,n_ir,ir(n_ir)
      REAL*8 :: v(3,n),scale
c      REAL   , INTENT(IN) :: v(3,n)

      INTEGER npts,iobj,i1,ilist,nlist,list(1000000)
      LOGICAL addobj
      REAL*8  a(3),b(3),p(3),pdist,t

      npts = n


      SELECTCASE(1) 
        CASE (1)

          IF (.FALSE.) THEN
            npts = 2
c            r  (1:2) = 0.5D0
c            phi(1:2) = 10.0D0
c            v(1,1:npts) =  r(1:npts) * DCOS(phi(1:npts)*D_DEGRAD)
c            v(3,1:npts) =  r(1:npts) * DSIN(phi(1:npts)*D_DEGRAD)
c            v(2,1) =  1.0D0
c            v(2,2) = -1.0D0

c            v(1,1) =  0.2D0
c            v(2,1) =  1.0D0
c            v(3,1) =  0.0D0
c            v(1,2) =  0.8D0
c            v(2,2) = -1.0D0
c            v(3,2) =  0.6D0

            v(1,1) =  0.5D0
            v(2,1) =  0.0D0
            v(3,1) =  0.0D0
            v(1,2) =  2.0D0
            v(2,2) =  0.0D0
            v(3,2) =  0.0D0

            WRITE(0,*) 'x:',v(1,1:npts)
            WRITE(0,*) 'y:',v(2,1:npts)
            WRITE(0,*) 'z:',v(3,1:npts)
          ENDIF

c          p(1) =  1.0D0
c          p(2) =  0.9999999D0
c          p(3) =  0.0D0
c          pdist = CalcPerp(v(1,1),v(1,2),p,t)

          CALL CalcDerivedQuantity(MODE_OBJ_CENTRE)

c...      Make this here for now, but it needs to be pre-computed:         
          nlist = 0
          DO iobj = 1, nobj
            IF (grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID) CYCLE

c            IF (obj(iobj)%index(IND_IR).NE.17) CYCLE
c            IF (obj(iobj)%index(IND_IR).NE.5) CYCLE

            addobj = .FALSE.  ! Need to completely redo this, just taking a range of items from precomputed lists for 
            DO i1 = 1, n_ir   ! each ring...
              IF (ir(i1).EQ.obj(iobj)%index(IND_IR)) addobj = .TRUE.
            ENDDO
            IF (addobj) THEN
              nlist = nlist + 1
              IF (nlist.EQ.1000000) STOP 'LIST: BUSTED!'
              list(nlist) = iobj
            ENDIF
          ENDDO


          DO ilist = 1, nlist
            iobj = list(ilist)

            IF (MOD(ilist,10000).EQ.0) WRITE(0,*) '  PROCESSING:',ilist

            IF (iobj.LT.10) WRITE(0,*) 'CEN:',obj_centroid(1:3,iobj)

            IF (obj(iobj)%segment(1).GT.0) CYCLE

            p(1:3) = obj_centroid(1:3,iobj)

            DO i1 = 1, npts-1
              a(1:3) = v(1:3,i1  )
              b(1:3) = v(1:3,i1+1)
c              a(1:3) = DBLE(v(1:3,i1  ))
c              b(1:3) = DBLE(v(1:3,i1+1))
              pdist = CalcPerp(a,b,p,t)
c              IF (pdist.NE.-1.0D0) THEN
c                WRITE(0,*) 'PDIST=',i1,iobj,pdist
c                WRITE(0,*) '     =',v(1:3,i1 )
c                WRITE(0,*) '     =',v(1:3,i1+1)
c              ENDIF
              IF (pdist.GE.0.0D0.AND.pdist.LE.scale) THEN
                obj(iobj)%segment(1) = 1  ! *** HACK ***
c                WRITE(0,*) '   BUSTED:',iobj
c                WRITE(0,*) '         :',obj_centroid(1:3,iobj)
                EXIT
c              ELSE
c                obj(iobj)%segment(1) = 0
              ENDIF
            ENDDO
          ENDDO

          CALL ClearDerivedQuantity(MODE_OBJ_CENTRE)

      ENDSELECT

      RETURN
 99   STOP
      END
