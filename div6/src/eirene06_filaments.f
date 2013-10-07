c     -*-Fortran-*-
c
c
c ======================================================================
c
c
c
c
      SUBROUTINE AssignFilamentPlasma
      USE mod_geometry
      USE mod_filament_params
      USE mod_filament
      USE mod_eirene06_locals
      USE mod_eirene06  ! *** TEMP ***
      USE mod_options
      IMPLICIT none

      REAL FindSeparatrixRadius

      INTEGER nplasma,iobj,count,iplasma,icell,ivtx,ifilament,
     .        nplasma_start,iplasma1,iplasma2

      INTEGER mode,nir,nmid,i1,i2,chop
      LOGICAL found_data
      REAL    x,y,z,rsep,r,
     .        t,frac,ne,te

      INTEGER n,index(10000),ring
      REAL    fraction(10000)
      REAL*8  len1,len2,v(3,10000),scale


      INTEGER, ALLOCATABLE :: tmp_segment(:)
      REAL   , ALLOCATABLE :: tmp_plasma(:,:)      


c      TYPE(type_object) newobj
c      IF (.NOT.ALLOCATED(plasma)) THEN
c        i1 = AddObject(newobj)
c        ALLOCATE(plasma(20,nobj))
c      ENDIF


c...  Store flag that marks the locally resolved tetrahedrons:
      ALLOCATE(tmp_segment(nobj))
      tmp_segment(1:nobj) = 0
c      tmp_segment(1:nobj) = obj(1:nobj)%segment(1)

c...  Find current size of PLASMA array:
      nplasma = 0
      DO iobj = 1, nobj
        nplasma = MAX(nplasma,obj(iobj)%index(IND_PLASMA))
      ENDDO   
      nplasma_start = nplasma

      CALL CalcDerivedQuantity(MODE_OBJ_CENTRE)
      CALL CalcDerivedQuantity(MODE_OBJ_TUBE)      
      CALL CalcDerivedQuantity(MODE_OBJ_DISTANCE)

      DO ifilament = 1, nfilament
        IF (filament(ifilament)%status.EQ.0) CYCLE

        obj_tube = 0
        obj_distance = 1.0D+20

        obj(1:nobj)%segment(1) = 0  

c...    Identify objects associated with this filament:
        mode = 1              ! Just compare perpendicular distance between 
        scale = opt_fil%scale(3) ! 0.005D0 ! 0.01D0        ! line segment and tetrahedron centroid 
                              ! For MODE=1, distance threshold for selection
        nir = 1 ! filament(ifilament)%ir_space(0,icell)
        DO icell = 1, filament(ifilament)%ncell
          ivtx = icell  ! *** HACK ***
          x = SNGL(filament(ifilament)%vtx(1,ivtx))
          y = 0.0
          z = SNGL(filament(ifilament)%vtx(3,ivtx))
          len1 = filament(ifilament)%lcell(1,ivtx) 
          len2 = filament(ifilament)%lcell(2,ivtx) 
          chop = 4
          IF (opt_fil%clip.EQ.1) chop = 6
          IF (opt_fil%clip.EQ.2) THEN
            len1 = 1.0D+20
            len2 = 0.0D0
            rsep = FindSeparatrixRadius(1)
            r = SQRT(x**2 + z**2)
            IF (r.LE.rsep) chop = 6
c            IF (chop.EQ.6) WRITE(0,*) '====WORKING:',r,rsep,
c     .         filament(ifilament)%ir_space(1,icell)
          ENDIF
          CALL TraceFieldLine_DIVIMP(x,y,z,2,chop,len1,len2,0.0D0,
     .                               n,v,index,fraction,ring,10000)  ! *** HACK *** (the 10000)
          CALL SelectTetrahedrons(n,v(1:3,1:n),mode,scale,
     .           ifilament,icell)
c     .           nir,filament(ifilament)%ir_space(1:nir,icell))
        ENDDO  ! TUBE


c...    Calculate the plasma quantities (as a function of time:
        t = SNGL(filament(ifilament)%t_last)

c...    Density:
        SELECTCASE (filament(ifilament)%ne_opt)
          CASE(2)
            n = (NINT(filament(ifilament)%ne_param(1)))
            ne = 0.0
            DO i1 = 2, n
              found_data = .FALSE.
              IF ((t+1.D-07.GE.filament(ifilament)%ne_param(i1  )).AND.
     .            (t-1.D-07.LE.filament(ifilament)%ne_param(i1+1))) THEN
                found_data = .TRUE.
                frac = (t - filament(ifilament)%ne_param(i1)) / 
     .                 (filament(ifilament)%ne_param(i1+1) - 
     .                  filament(ifilament)%ne_param(i1  ))
                i2 = i1 + n
                ne = (1.0D0-frac)*filament(ifilament)%ne_param(i2  ) +
     .                      frac *filament(ifilament)%ne_param(i2+1)
              ENDIF
            ENDDO
          CASE DEFAULT
            STOP 'PROBLEM XXX1'
        ENDSELECT
        WRITE(0,*) 't,ne:',ifilament,NINT(t/1.0D-06),ne,frac

c...    Temperature:
        SELECTCASE (filament(ifilament)%te_opt)
          CASE(2)
            n = (NINT(filament(ifilament)%te_param(1)))
            te = 0.0
            DO i1 = 2, n
              found_data = .FALSE.
              IF ((t+1.D-07.GE.filament(ifilament)%te_param(i1  )).AND.
     .            (t-1.D-07.LE.filament(ifilament)%te_param(i1+1))) THEN
                found_data = .TRUE.
                frac = (t - filament(ifilament)%te_param(i1)) / 
     .                 (filament(ifilament)%te_param(i1+1) - 
     .                  filament(ifilament)%te_param(i1  ))
                i2 = i1 + n
                te = (1.0D0-frac)*filament(ifilament)%te_param(i2  ) +
     .                      frac *filament(ifilament)%te_param(i2+1)
              ENDIF
            ENDDO
          CASE DEFAULT
            STOP 'PROBLEM XXX1'
        ENDSELECT
        WRITE(0,*) 't,te:',ifilament,t/1.0D-06,te,frac

c...    Increase size of PLASMA array:
        ALLOCATE(tmp_plasma(20,nplasma))  ! Needs to be consistent with PLASMA allocation in ProcessTetrahedrons
        tmp_plasma = plasma
        DEALLOCATE(plasma)

c...    Find the number of associated cells:
        count = 0
        DO iobj = 1, nobj
          IF (obj(iobj)%segment(1).EQ.1.AND.
     .        obj(iobj)%index(IND_PLASMA).LE.nplasma_start)     
     .      count = count + 1
        ENDDO
        WRITE(0,*) 'nplasma,count:',nplasma,count
c        count = 1  ! *** HACK *** 

        ALLOCATE(plasma(20,nplasma+count))  ! Needs to be consistent with PLASMA allocation in ProcessTetrahedrons
        DO iplasma = 1, nplasma
          plasma(1:20,iplasma) = tmp_plasma(1:20,iplasma)
        ENDDO
        DEALLOCATE(tmp_plasma)

        count = 0
        DO iobj = 1, nobj                         ! Fast if I kept a list of objects from the above loop...
          IF (obj(iobj)%segment(1).EQ.0) CYCLE  

          iplasma1 = obj(iobj)%index(IND_PLASMA)
          IF (iplasma1.LE.nplasma_start) THEN
            count = count + 1
            iplasma2 = nplasma + count
          ELSE
            iplasma2 = iplasma1
          ENDIF
          
          frac = ne / (ne + plasma(3,iplasma1))

c          IF (iplasma2.EQ.iplasma1)
c     .      WRITE(0,*) '-OLD>',plasma(1:3,iplasma2)
c          IF (iplasma2.NE.iplasma1)
c     .      WRITE(0,*) '-NEW>',plasma(1:3,iplasma2)

          plasma(1:2,iplasma2) = (1.0-frac) * plasma(1:2,iplasma1) + 
     .                                frac  * te
          plasma(3  ,iplasma2) = plasma(3,iplasma1) + ne
          plasma(17:19,iplasma2) = plasma(17:19,iplasma1)  ! *** HACK *** For local scaling otetrahedron target flux...
c          plasma(1,nplasma) = te       ! Te (eV)
c          plasma(2,nplasma) = te       ! Ti (eV)
c          plasma(3,nplasma) = ne       ! ni (eV) (ne=ni assumed at present)

c          IF (iplasma2.EQ.iplasma1)
c     .      WRITE(0,*) '     ',plasma(1:3,iplasma1),frac
c          IF (iplasma2.EQ.iplasma1)
c     .      WRITE(0,*) '     ',plasma(1:3,iplasma2)
c          IF (iplasma2.NE.iplasma1)
c     .      WRITE(0,*) '     ',plasma(1:3,iplasma1),frac
c          IF (iplasma2.NE.iplasma1)
c     .      WRITE(0,*) '     ',plasma(1:3,iplasma2),
c     .  obj(iobj)%index(IND_IK),obj(iobj)%index(IND_IR)


          obj(iobj)%index(IND_PLASMA) = iplasma2

          tmp_segment(iobj) = 1
        ENDDO           

        nplasma = nplasma + count
c        WRITE(0,*) 'nplasma,count:',nplasma,count

      ENDDO  ! FILAMENT

      CALL ClearDerivedQuantity(MODE_OBJ_CENTRE)
      CALL ClearDerivedQuantity(MODE_OBJ_TUBE)
      CALL ClearDerivedQuantity(MODE_OBJ_CENTRE)

      obj(1:nobj)%segment(1) = tmp_segment(1:nobj)
      DEALLOCATE(tmp_segment)

c        STOP 'dggsdgsd'



      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
c *** Need to look over zones maybe since each fully split tetrahedron
c     goes into 24 tetrahedrons, which may blow the gasket *** 
c
c
      SUBROUTINE ResolveFilament(itet)
      USE mod_geometry
      USE mod_filament_params
      USE mod_filament
      USE mod_eirene06_locals
      USE mod_eirene06  ! *** TEMP ***
      USE mod_options
      IMPLICIT none

      INTEGER, INTENT(IN) :: itet

c      REAL FindSeparatrixRadius
      REAL*8 CalcPerp

      INTEGER iobj,iobj1,tmp_nobj,ifilament,ivtx,icell,iside,mode,iloop

      INTEGER nir,nmid,chop
      REAL    x,y,z

      INTEGER n,index(10000),ring,opt_resolve
      REAL    fraction(10000)
      REAL*8  len1,len2,v(3,10000),scale,a(3),b(3),p(3),t,pdist

      IF (itet.EQ.-1) THEN
        opt_resolve = 3
      ELSE
        opt_resolve = opt_eir%tet_mode(itet)
      ENDIF

      obj(1:nobj)%segment(1) = 0

c      RETURN

      SELECTCASE (opt_resolve)
c      SELECTCASE (4)
c       ----------------------------------------------------------------
        CASE(-1)
c       ----------------------------------------------------------------
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
c       ----------------------------------------------------------------
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
c       ----------------------------------------------------------------
        CASE (2)
          CALL TraceFieldLine_DIVIMP(0.0,0.0,0.0,1,1,1.0E+20,1.0E+20,
     .                              0.0D0,n,v,index,fraction,ring,10000)
          CALL SelectTetrahedrons(n,v(1:3,1:n))
          tmp_nobj = nobj
          DO iobj = 1, tmp_nobj
            IF (obj(iobj)%segment(1).GT.0) CALL DivideTetrahedron(iobj)
          ENDDO
          CALL CleanObjectArray
c       ----------------------------------------------------------------
        CASE (3)

          DO iloop = 1, 2
            obj(1:nobj)%segment(1) = 0

            mode = 1              ! Just compare perpendicular distance between 
            IF (iloop.EQ.1) THEN  ! line segment and tetrahedron centroid 
              scale = opt_fil%scale(1) ! 0.01D0   ! 0.05D0      ! For MODE=1, distance threshold for selection
            ELSE
              scale = opt_fil%scale(2) ! 0.005D0  ! 0.02D0
            ENDIF
            nir = 1 ! filament(ifilament)%ir_space(0,icell)

            CALL CalcDerivedQuantity(MODE_OBJ_CENTRE)
            DO ifilament = 1, nfilament
              IF (filament(ifilament)%status.EQ.0) CYCLE

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
c                WRITE(0,*) '==FILAMENT:',ifilament,ivtx,x,z,
c     .            nir,filament(ifilament)%ir_space(1:nir,icell)
                chop = 4
                IF (opt_fil%clip.EQ.1) chop = 6
                IF (opt_fil%clip.EQ.2) THEN
                  len1 = 1.0D+20
                  len2 = 0.0D0
                ENDIF
                IF (opt_fil%clip.EQ.3) THEN
                  len1 = 1.0D+20
                  len2 = 1.0D+20
                ENDIF
c         WRITE(0,*) '============>LENGTH:',len1,len2,chop
                CALL TraceFieldLine_DIVIMP(x,y,z,2,chop,len1,len2,0.0D0,
     .                                   n,v,index,fraction,ring,10000)  ! *** HACK *** (the 10000)
                CALL SelectTetrahedrons(n,v(1:3,1:n),mode,scale,
     .                 ifilament,icell)
c     .                 nir,filament(ifilament)%ir_space(1:nir,icell))
              ENDDO  ! TUBE
            ENDDO  ! FILAMENT
            CALL ClearDerivedQuantity(MODE_OBJ_CENTRE)

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
                IF (obj(iobj)%segment(1).GT.0)         ! *** BUG *** do do this only on the first pass?  -SL, 29/03/2010
     .            CALL DivideTetrahedron(iobj)
              ENDDO
            ENDIF

            CALL CleanObjectArray
          ENDDO  ! LOOP
c       ----------------------------------------------------------------
        CASE (4)  ! Resolve if near a specified line
          a(1) = opt_eir%tet_x1(itet)
          a(2) = opt_eir%tet_y1(itet)
          a(3) = opt_eir%tet_z1(itet)
          b(1) = opt_eir%tet_x2(itet)
          b(2) = opt_eir%tet_y2(itet)
          b(3) = opt_eir%tet_z2(itet)

c          write(0,*) 'resolving... a=',a(1:3)
c          write(0,*) 'resolving... b=',b(1:3)

          DO iloop = 1, NINT(opt_eir%tet_param1(itet))
c          DO iloop = 1, 1

            obj(1:nobj)%segment(1) = 0
            IF     (iloop.LE.1) THEN  ! line segment and tetrahedron centroid 
              scale = opt_eir%tet_param2(itet)  ! 0.050D0 ! 0.100D0   
c           ELSEIF (iloop.LE.2) THEN  ! line segment and tetrahedron centroid 
c             scale = 0.150D0   
            ELSE
              scale = scale / 2.0D0 
c              scale = 0.020D0 
            ENDIF

            WRITE(0,*) 'resolving... pass=',iloop,scale


            CALL CalcDerivedQuantity(MODE_OBJ_CENTRE)

            DO iobj = 1, nobj                       
              p(1:3) = obj_centroid(1:3,iobj)
              pdist = CalcPerp(a,b,p,t)
c              IF (pdist.NE.-1.0D0) THEN
c                WRITE(0,*) 'PDIST=',iobj,pdist
c                WRITE(0,*) '     =',a(1:3)
c                WRITE(0,*) '     =',b(1:3)
c              ENDIF
              IF (pdist.GE.0.0D0.AND.pdist.LE.scale) THEN
                obj(iobj)%segment(1) = 1  ! *** HACK ***
              ENDIF
            ENDDO
            CALL ClearDerivedQuantity(MODE_OBJ_CENTRE)

            IF (iloop.EQ.1) THEN
c...          Resolve the nighbouring tetrahedrons as well to make sure that the
c             region of interest is properly resolved:
              DO iobj = 1, nobj                       
                IF (obj(iobj)%segment(1).NE.1) CYCLE  
                DO iside = 1, obj(iobj)%nside         
                  iobj1 = obj(iobj)%omap(iside)
                  IF (iobj1.LE.0) CYCLE
                  IF (obj(iobj1)%segment(1).EQ.0) 
     .              obj(iobj1)%segment(1) = 2  ! *** HACK ***
                ENDDO
              ENDDO           
            ENDIF

c...        Increase spatial resolution:
            tmp_nobj = nobj
            DO iobj = 1, tmp_nobj
c              WRITE(0,*) 'what',obj(iobj)%segment(1)
c              IF (obj(iobj)%segment(1).GT.0)
c     .    WRITE(0,*) ' DIVIDING TETRAHEDRON!',iobj 
              IF (obj(iobj)%segment(1).GT.0) 
     .          CALL DivideTetrahedron(iobj)
            ENDDO

            CALL CleanObjectArray

            write(0,*) 'old/new',tmp_nobj,nobj,nobj-tmp_nobj

          ENDDO  ! LOOP
c       ----------------------------------------------------------------
      ENDSELECT

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
      SUBROUTINE SelectTetrahedrons(n,v,mode,scale,ifilament,icell)
c      SUBROUTINE SelectTetrahedrons(n,v,mode,scale,n_ir,ir)
      USE mod_geometry
      USE mod_filament
      USE mod_eirene06_locals
      IMPLICIT none

      REAL*8 CalcPerp

      INTEGER, INTENT(IN) :: n,mode,ifilament,icell
      REAL*8 :: v(3,n),scale
c      REAL   , INTENT(IN) :: v(3,n)

      INTEGER npts,iobj,i1,ilist,nlist,list(1000000),n_ir,ir(10)
      LOGICAL addobj
      REAL*8  a(3),b(3),p(3),pdist,t,distance(n),delta,length,len1

      npts = n

      n_ir = 1
      ir(1:n_ir) = filament(ifilament)%ir_space(1:n_ir,icell)

      IF (.NOT.ALLOCATED(obj_centroid)) 
     .  CALL ER('SelectTetrahedrons','Need centroids',*99)


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

c...      Calculate the distance along the field line for each
c         segment, using the midplane as 0.0:
          len1 = filament(ifilament)%lcell(1,icell)
          DO i1 = 1, npts-1
            delta =  DSQRT((v(1,i1+1)-v(1,i1))**2 + 
     .                     (v(2,i1+1)-v(2,i1))**2 + 
     .                     (v(3,i1+1)-v(3,i1))**2)
            length = length + 0.5D0 * delta
            distance(i1) = -len1 + length
            length = length + 0.5D0 * delta       
          ENDDO
c          WRITE(0,*) '=====> LENGTH:',length,
c     .      len1 + filament(ifilament)%lcell(1,icell)

          DO ilist = 1, nlist
            iobj = list(ilist)

c            IF (MOD(ilist,1000000).EQ.0) 
c     .        WRITE(0,*) '  PROCESSING:',ilist

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

c                IF (pdist.LT.obj_distance(1,iobj)) THEN  ! *** LEFT OFF ***
c                  obj_tube = icell
c                  obj_distance(1,iobj) = pdist
c                  obj_distance(2,iobj) = distance(i1)
c                ENDIF

                obj(iobj)%segment(1) = 1  ! *** HACK ***
c                WRITE(0,*) '   BUSTED:',iobj
c                WRITE(0,*) '         :',obj_centroid(1:3,iobj)
                EXIT
c              ELSE
c                obj(iobj)%segment(1) = 0
              ENDIF
            ENDDO
          ENDDO


      ENDSELECT

      RETURN
 99   STOP
      END
