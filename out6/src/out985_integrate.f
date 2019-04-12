c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE GenerateReflectionChords(mode,imodel,weight,vv,nv,
     .                                    nref,rv,rw,sw,reflvl)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

c...  Input:
      INTEGER mode,imodel,nref,i,j,i1,reflvl
      REAL*8 weight,vv(3),nv(3),rv(3,nref),rw(nref),sw(nref),ra(nref),
     .       mat1(3,3),mat2(3,3),mat3(3,3),
     .       theta,phi,dtheta,dphi,solid_angle,solid_total,
     .       cutoff,dtheta1,dphi1,scale_total,scale_gross

      INTEGER iphi,nphi,reflection_model,iref,nexp
      REAL*8 ndotv,ndotr,v1(3),angle,ddum1,scale_factor

      rw = 0.0D0

      reflection_model = opt%ref_model(imodel)

      scale_gross = weight

      SELECTCASE (opt%ref_wlgth(imodel))
        CASE(0) 
          scale_gross = scale_gross * DBLE(opt%ref_k(imodel)**reflvl)       ! LEVEL=2: need to square this... 
        CASE DEFAULT
          CALL User_SurfaceReflectivity(imodel,scale_gross)
      ENDSELECT
c      WRITE(0,*) 'scale_gross=',scale_gross

      IF (output) WRITE(0,*) 'REFLECTION MODEL:',reflection_model

      IF (.TRUE.) THEN

        SELECTCASE (reflection_model)
          CASE (1)
            dtheta1 =   DBLE(opt%ref_dtheta(imodel))
            dphi1   =   DBLE(opt%ref_dphi  (imodel))
            nexp    =   opt%ref_n(imodel)
            cutoff  =   DBLE(opt%ref_cutoff(imodel))
          CASE(2)
            dtheta1 =   DBLE(opt%ref_dtheta(imodel))
            dphi1   =   DBLE(opt%ref_dphi  (imodel))
            nexp    =   1
            cutoff  =  -1.0D+10
          CASE DEFAULT
            nref = 0
            RETURN
        ENDSELECT

        IF (opt%ref_ophi(imodel).EQ.0) THEN  ! Turn off rays outside poloidal plane
          dphi1 = 180.0
        ENDIF

        dtheta = dtheta1
        dphi   = dphi1

c...    Primary reflection chord:
        nref = 1        
c        solid_angle = dtheta * D_DEGRAD
        IF (opt%ref_ophi(imodel).EQ.0) THEN
          solid_angle = dtheta * D_DEGRAD
        ELSE
          solid_angle = DCOS(0.0D0) - 
     .                  DCOS(0.5D0*dtheta*D_DEGRAD)
        ENDIF
        scale_factor = 1.0D0
        rv(1,1) = 0.0D0
        rv(2,1) = 0.0D0
        rv(3,1) = 1.0D0
        rw(1) = scale_factor * solid_angle
        sw(1) = scale_factor
        ra(1) = solid_angle

        IF (output) WRITE(0,*) '  DATA:',1,solid_angle

c...    Loop to generate associated reflection chords:
        theta = 0.0D0

        DO WHILE (.TRUE.) 

          theta = theta + 0.5D0 * dtheta

          SELECTCASE (opt%ref_otheta(imodel))
            CASE (0)
              EXIT
            CASE (1)
              dtheta = dtheta1
            CASE DEFAULT
              STOP 'NOT TO BE HERE'
          ENDSELECT

          theta = theta + 0.5D0 * dtheta

c          scale_factor = 1.0D0
          scale_factor = DCOS(theta*D_DEGRAD)**nexp

          IF (scale_factor.LT.cutoff.OR.theta+dtheta.GT.90.0D0) EXIT       ! Exit conditions   

c          IF (output) WRITE(0,*) '  DATA:',nref,scakele_factor,theta


          dphi = dphi1
          nphi = NINT(360.0D0 / dphi)

          IF (opt%ref_ophi(imodel).EQ.0) THEN
            solid_angle = 2.0D0 * dtheta * D_DEGRAD
          ELSE
            solid_angle = DCOS((theta-0.5D0*dtheta)*D_DEGRAD) - 
     .                    DCOS((theta+0.5D0*dtheta)*D_DEGRAD)
          ENDIF

          IF (nref+nphi.GT.10000000) 
     .      CALL ER('Reflections','Exceeded limit',*99)


          solid_angle = solid_angle / DBLE(REAL(nphi))

c          WRITE(0,'(A,I6,5F12.5)') 
c     .      '  DATA:',nref,scale_factor,theta,
c     .                solid_angle/ DBLE(REAL(nphi)),
c     .           scale_factor,solid_angle

          nref = nref + 1
          rv(1:3,nref) = rv(1:3,1)
          rw(nref) = scale_factor * solid_angle
          sw(nref) = scale_factor
          ra(nref) = solid_angle
c...
          CALL Calc_Transform2(mat1,0.0D0,1,0)
          CALL Calc_Transform2(mat1,theta*D_DEGRAD,1,1)
          CALL Transform_Vect (mat1,rv(1,nref))
c...      
          CALL Calc_Transform2(mat1,0.0D0,3,0)
          CALL Calc_Transform2(mat1,dphi*D_DEGRAD,3,1)

c...      *** Ideally, want another rotation here I think, to align this first vector
c             with the incident view, so that there is the option to weight
c             the reflecitons in the plane of the normal+incident vector... ***

          DO iphi = 2, nphi       ! *** problem if dphi doesn't go into 360.0 evenly... ***
            nref = nref + 1
            rv(1:3,nref) = rv(1:3,nref-1)
c...
            CALL Transform_Vect(mat1,rv(1,nref))
c...        Calculate solid angle and set weight:
            rw(nref) = scale_factor * solid_angle
            sw(nref) = scale_factor
            ra(nref) = solid_angle

            IF (output) THEN
              WRITE(0,*) 'RV:',nref,rv(1:3,nref)
            ENDIF   
          ENDDO
        ENDDO
      ENDIF

      SELECTCASE (reflection_model)
        CASE (1)  ! Specular
c...      Calculate specular reflection vector:
          ndotv = nv(1) * vv(1) + nv(2) * vv(2) + nv(3) * vv(3)
          v1(1:3) = 2.0D0 * ndotv * nv(1:3) - vv(1:3)
          rv(1:3,1) = v1(1:3)
        CASE (2)  ! Diffuse
          v1(1:3) = nv(1:3)
        CASE DEFAULT
           STOP 'STOP: SOMETHING WRONG'
      ENDSELECT


      IF (output) THEN
        WRITE(0,*) 'VV:',vv(1:3)
        WRITE(0,*) 'NV:',nv(1:3)         
        WRITE(0,*) 'V1:',v1(1:3)           
        WRITE(0,*) 'RV1:',rv(1:3,1)
        WRITE(0,*) 'ANGLE:',angle*180.0/3.1415
      ENDIF

! this should be for all reflection models, but it's not working, I don't think
      ! so turning it off for now 
      IF (reflection_model.EQ.2) THEN
c...    Rotate array of reflection vectors:

        stop 'not sure this works...' ! 05/01/2019
         
c       IF (v1(1).LT.0.0D0) THEN
c         angle = +1.0D0 * DATAN(v1(2) / (v1(1) + 1.0D-10))
c       ELSE
c         angle = -1.0D0 * DATAN(v1(2) / (v1(1) + 1.0D-10))
c       ENDIF

        angle = -1.0D0 * DASIN(v1(2))

c...    Rotate about x-axis (tilt):                       !... better to do swing before tilt? 
        CALL Calc_Transform2(mat1,0.0D0,1,0)
        CALL Calc_Transform2(mat1,angle,1,1)
c        CALL Transform_Vect(mat1,rv(1,1))

        IF (v1(1).LT.0.0D0) THEN
          angle = +1.0D0 * DACOS(v1(3))
        ELSE
          angle = -1.0D0 * DACOS(v1(3))
        ENDIF

        IF (output) WRITE(0,*) 'ANGLE:',angle*180.0/3.1415

c...    Rotate about y-axis (....):                       !... better to do swing before tilt? 
        CALL Calc_Transform2(mat2,0.0D0,2,0)
        CALL Calc_Transform2(mat2,angle,2,1)
c        CALL Transform_Vect(mat2,rv(1,1))

c...    Combine rotations:
        mat3 = 0.0D0
        DO i = 1, 3
          DO j = 1, 3
            mat3(i,j) = mat2(i,1) * mat1(1,j) + 
     .                  mat2(i,2) * mat1(2,j) + 
     .                  mat2(i,3) * mat1(3,j) 
          ENDDO
        ENDDO
        DO i1 = 1, nref
          CALL Transform_Vect(mat3,rv(1,i1))
        ENDDO


c        rv(1:3,1) = v1(1:3)
     
        IF (output) WRITE(0,*) 'RV:',rv(1:3,1)

      ENDIF


      IF (.TRUE.) THEN  ! *** THIS COULD BE DONE ABOVE, MAYBE FASTER, IF THE ORIENTATION OF r RELATIVE TO v/n WERE KNOWN
c...    Delete chords that are greater than 90 degrees from the normal, an issue
c       for specular reflections:

        DO iref = 1, nref
          ndotr = nv(1)*rv(1,iref) + nv(2)*rv(2,iref) + nv(3)*rv(3,iref)
          angle = DACOS(ndotr)
          IF (angle/D_DEGRAD.GT.90.0D0) rw(iref) = -999.0D0  ! Tag for deletion
        ENDDO

c...    Clean up:
        DO iref = nref, 1, -1
          IF (rw(iref).EQ.-999.0D0) THEN
            DO i1 = iref, nref-1
              rv(1:3,i1) = rv(1:3,i1+1)
              rw(    i1) = rw(    i1+1)
              sw(    i1) = sw(    i1+1)
            ENDDO
            nref = nref - 1
          ENDIF
        ENDDO

      ENDIF


      IF (.TRUE.) THEN

c...    Normalize total solid angle:
        solid_total = 0.0D0
        DO iref = 1, nref
          solid_total = solid_total + ra(iref)
c          WRITE(0,*) 'SA:',iref,ra(iref),solid_total
        ENDDO
        rw(1:nref) = rw(1:nref) / solid_total

c        DO iref = 1, nref
c          scale_total = scale_total + rw(iref)
c        ENDDO
c        rw(1:nref) = rw(1:nref) / scale_total

        IF (output) THEN
          ddum1 = 0.0D0
          DO i1 = 1, nref
            ddum1 = ddum1 + rw(i1)
          ENDDO
          WRITE(0,*) '  WEIGHT SUM:',nref,ddum1,solid_total
        ENDIF


        SELECTCASE (1)
          CASE(1)
            scale_factor = scale_gross
          CASE(2)
            scale_factor = scale_gross / rw(1)   ! Fixed fraction of central reflection chord is reflected
          CASE DEFAULT
        ENDSELECT

        rw(1:nref) = rw(1:nref) * scale_factor 
      ENDIF

c      WRITE(6,*) 'SCALE REF:',solid_total,scale_factor
c      WRITE(6,*) '         :',rw(1:nref)
c      WRITE(6,*) '         :',ra(1:nref)

        IF (output) WRITE(0,*) '  WEIGHT 1:',rw(1)
        IF (output) WRITE(0,*) '  WEIGHT  :',rw(2:nref)


c      WRITE(0,*) 'NREF:',nref

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE AssignReflectionChord(chord,iobj,iside,isrf,
     .                                 refcnt,reflvl)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

c...  Input:
      INTEGER iobj,iside,isrf,refcnt,reflvl
      TYPE(type_view) :: chord

      INTEGER nref,imodel
      REAL*8  av(3),bv(3),nv(3),vv(3),tv(3),length,ndotv,theta

      REAL*8, ALLOCATABLE :: rv(:,:),rw(:),sw(:)

      TYPE(type_view) :: chord_hold,chord_primary
 
      SAVE

      IF (output) WRITE(0,*) 'REFLECTING:',iobj,iside,isrf

c      RETURN

c     ------------------------------------------------------------------
      IF (refcnt.EQ.1) THEN
c...    Initialize:    

c...    Store incident chord view+weight, normalize:
        IF (reflvl.EQ.1) chord_primary = chord

        chord_hold = chord
c        WRITE(0,*) 'INTEGRAL:',nchord,chord%integral(1)

        IF (output) THEN
          WRITE(0,'(A,3F10.3)') 'ref  chord%v1,2=',chord%v1
          WRITE(0,'(A,3F10.3)') '               =',chord%v2
        ENDIF
          
        vv(1:3) = chord%v1(1:3) - chord%v2(1:3)
c...    Normalize:
        length = DSQRT(vv(1)**2 + vv(2)**2 + vv(3)**2) 
        vv(1:3) = vv(1:3) / length

c...    Calcualte surface normal, normalize:
        IF     (obj(iobj)%gsur(iside).EQ.GT_TC) THEN
c....       
          IF (DABS(chord%v2(1)).GT.1.0D-10) THEN
            theta = DATAN(chord%v2(3) / chord%v2(1))
          ELSE
            theta = 90.0D0 * D_DEGRAD
          ENDIF

          IF (obj(iobj)%nside.NE.0) THEN

            av(1) = vtx(1,srf(isrf)%ivtx(1))*
     .              DABS(DCOS(theta)) * DSIGN(1.0D0,chord%v2(1))
            av(2) = vtx(2,srf(isrf)%ivtx(1))
            av(3) = vtx(1,srf(isrf)%ivtx(1))*
     .              DABS(DSIN(theta)) * DSIGN(1.0D0,chord%v2(3))

            tv(1) = vtx(1,srf(isrf)%ivtx(2))*
     .              DABS(DCOS(theta)) * DSIGN(1.0D0,chord%v2(1))
            tv(2) = vtx(2,srf(isrf)%ivtx(2))
            tv(3) = vtx(1,srf(isrf)%ivtx(2))*
     .              DABS(DSIN(theta)) * DSIGN(1.0D0,chord%v2(3))

          ELSE

            IF (obj(iobj)%ipts(1,iside).EQ.0.OR.
     .          obj(iobj)%ipts(1,iside).EQ.0) THEN
              STOP 'PROBLEM AAA2'
            ENDIF

            av(1) = obj(iobj)%v(1,obj(iobj)%ipts(1,iside)) *
     .              DABS(DCOS(theta)) * DSIGN(1.0D0,chord%v2(1))
            av(2) = obj(iobj)%v(2,obj(iobj)%ipts(1,iside))
            av(3) = obj(iobj)%v(1,obj(iobj)%ipts(1,iside)) *
     .              DABS(DSIN(theta)) * DSIGN(1.0D0,chord%v2(3))

            tv(1) = obj(iobj)%v(1,obj(iobj)%ipts(2,iside)) *
     .              DABS(DCOS(theta)) * DSIGN(1.0D0,chord%v2(1))
            tv(2) = obj(iobj)%v(2,obj(iobj)%ipts(2,iside))
            tv(3) = obj(iobj)%v(1,obj(iobj)%ipts(2,iside)) *
     .              DABS(DSIN(theta)) * DSIGN(1.0D0,chord%v2(3))

c            WRITE(0,*) '  THETA=',theta,iobj,iside
c            WRITE(0,*) '  AV   =',SNGL(av)
c            WRITE(0,*) '  TV   =',SNGL(tv)
c            WRITE(0,*) '  V2   =',SNGL(chord%v2(1:3))

c
c            IF (DABS(chord%v2(1)).GT.1.0D-10) THEN
c              theta = DATAN(chord%v2(3) / chord%v2(1))
c            ELSE
c              theta = 90.0D0 * D_DEGRAD
c            ENDIF
c            av(1) = obj(iobj)%v(1,obj(iobj)%ipts(1,iside)) *
c     .              DABS(DCOS(theta)) * DSIGN(1.0D0,chord%v2(1))
c            av(2) = obj(iobj)%v(2,obj(iobj)%ipts(1,iside))
c            av(3) = obj(iobj)%v(1,obj(iobj)%ipts(1,iside)) *
c     .              DABS(DSIN(theta)) * DSIGN(1.0D0,chord%v2(3))
c            av(1:3) = av(1:3) - chord%v2(1:3)
c
c            bv(1) =  chord%v2(3)      ! Cross-product of the surface intersection
c            bv(2) =  0.0D0            ! vector with the y-axis unit vector
c            bv(3) = -chord%v2(1) 
          ENDIF

c          WRITE(0,*) '  THETA=',theta,iobj,iside
c          WRITE(0,*) '  TV   =',SNGL(tv)
c          WRITE(0,*) '  AV   =',SNGL(av)

          av(1:3) = av(1:3) - tv(1:3)
c          av(1:3) = av(1:3) - chord%v2(1:3)
c          WRITE(0,*) '  AV   =',SNGL(av)

          bv(1) =  chord%v2(3)      ! Cross-product of the surface intersection
          bv(2) =  0.0D0            ! vector with the y-axis unit vector
          bv(3) = -chord%v2(1) 
c          WRITE(0,*) '  BV   =',SNGL(bv)

c          WRITE(0,*) '  AV   =',SNGL(av)
c          WRITE(0,*) '  BV   =',SNGL(bv)

        ELSEIF (obj(iobj)%gsur(iside).EQ.GT_TD) THEN
c...      Floating 3D cartesian surface (most general toroidally 
c         discretized representation):

          IF (output) WRITE(0,*) 'GT_TD detected' 
           
c...      Store this?  Certainly...
          IF (obj(iobj)%nside.NE.0) THEN
            IF (srf(isrf)%ivtx(1).EQ.0.OR.
     .          srf(isrf)%ivtx(1).EQ.0) THEN
              STOP 'PROBLEM BBB2'
            ENDIF
            av(1:3) = vtx(1:3,srf(isrf)%ivtx(1)) - 
     .                vtx(1:3,srf(isrf)%ivtx(2))
            bv(1:3) = vtx(1:3,srf(isrf)%ivtx(3)) - 
     .                vtx(1:3,srf(isrf)%ivtx(2))
          ELSE
            IF (obj(iobj)%ipts(1,iside).EQ.0.OR.
     .          obj(iobj)%ipts(1,iside).EQ.0) THEN
              STOP 'PROBLEM AAA2'
            ENDIF
            av(1:3) = obj(iobj)%v(1:3,obj(iobj)%ipts(1,iside)) - 
     .                obj(iobj)%v(1:3,obj(iobj)%ipts(2,iside))
            bv(1:3) = obj(iobj)%v(1:3,obj(iobj)%ipts(3,iside)) - 
     .                obj(iobj)%v(1:3,obj(iobj)%ipts(2,iside))
          ENDIF
        ENDIF

c...    Surface normal:
        nv(1) =  av(2) * bv(3) - av(3) * bv(2)
        nv(2) = -av(1) * bv(3) + av(3) * bv(1)
        nv(3) =  av(1) * bv(2) - av(2) * bv(1) 
c...    Normalize:
        length = DSQRT(nv(1)**2 + nv(2)**2 + nv(3)**2) 
        nv(1:3) = nv(1:3) / length

c...    Make sure the correct normal is in use, relative to the 
c       viewing chord:
        theta = DACOS(vv(1)*nv(1) + vv(2)*nv(2) + vv(3)*nv(3)) /
     .          D_DEGRAD
        IF (output) THEN
          WRITE(0,*) 'THETA:',theta
          WRITE(0,*) 'vv:',REAL(vv) 
          WRITE(0,*) 'nv:',REAL(nv) 
        ENDIF
c          WRITE(0,*) 'THETA:',theta
c          WRITE(0,*) 'vv:',REAL(vv) 
c          WRITE(0,*) 'nv:',REAL(nv) 
        IF (theta.GT.90.0D0) nv(1:3) = -nv(1:3)
c          WRITE(0,*) 'nv:',REAL(nv) 
         
c...    Check model, calculate specular/diffuse reflection chord array:

c...    Call GenerateReflectionChords once to count how many...
c...    Call GenerateReflectionChords then again to assign...

        nref =10000000
        IF (ALLOCATED(rv)) THEN
          WRITE(0,*) '  UNEXPECTED RESET OF rv AND zv ARRAYS'
          DEALLOCATE(rv)
          DEALLOCATE(rw)
          DEALLOCATE(sw)
        ENDIF 

        ALLOCATE(rv(3,nref))
        ALLOCATE(rw(nref))
        ALLOCATE(sw(nref))

        imodel = obj(iobj)%reflec(iside)
        IF (imodel.GT.opt%ref_num) 
     .    CALL ER('GenerateReflectionChord','Reflection model '//     
     .            'appears invalid',*99)  

        CALL GenerateReflectionChords(0,imodel,chord%weight,vv,nv,
     .                                nref,rv,rw,sw,reflvl)

        IF (output) THEN
          WRITE(0,*) 'SURFACE:',iobj,iside,isrf
          WRITE(0,*) 'VIEW:   ',vv(1:3)
          WRITE(0,*) 'NORMAL: ',nv(1:3)
          WRITE(0,*) 'REF1:   ',rv(1:3,1)          
c          WRITE(0,*) 'REFLEC(1): ',reflec(1)%v2(1:3)
        ENDIF

        IF (output) WRITE(0,*) 'FINGERS CROSSED'

      ENDIF  ! End of initialization
c     ------------------------------------------------------------------

      IF (nref.EQ.0.OR.refcnt.GT.nref) THEN
        IF (output) WRITE(0,*) 'TURNING OFF REFECTION'
        refcnt = 0
        IF (ALLOCATED(rv)) DEALLOCATE(rv)
        IF (ALLOCATED(rw)) DEALLOCATE(rw)
        IF (ALLOCATED(sw)) DEALLOCATE(sw)
c...    Restore:
        chord%v1(1:3) = chord_primary%v1(1:3)  
        chord%v2(1:3) = chord_primary%v2(1:3) 
c        chord%v1(1:3) = chord_hold%v1(1:3)  
c        chord%v2(1:3) = chord_hold%v2(1:3) 
c        WRITE(0,*) 'INTEGRAL:',nchord,chord%integral(1)
        RETURN
      ENDIF

c       refcnt = 0

      chord%weight      = rw(refcnt)
      chord%weight_save = sw(refcnt)
      chord%v1(1:3) = chord_hold%v2(1:3)
      chord%v2(1:3) = chord_hold%v2(1:3) + 100.0D0 * rv(1:3,refcnt)

c      WRITE(0,*) '---------chord weight',chord%weight
c      WRITE(6,*) '---------chord weight',chord%weight

c      WRITE(0,*) 'CHORD LENGTH: ',
c     .  DSQRT((chord%v2(1)-chord%v1(1))**2 + 
c     .        (chord%v2(2)-chord%v1(2))**2 + 
c     .        (chord%v2(3)-chord%v1(3))**2)

      IF (output) WRITE(0,*) 'DONE MAKING REFLECTIONS'

      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE IntegrateChord(chord_primary,status)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

      TYPE(type_view)       :: chord_primary
      INTEGER status                    ! Flag nobj so it can't change?

      REAL*8     RTOL
      PARAMETER (RTOL=1.0D-06)

      INTEGER i1,i2,i3,iobj,isid,isrf,ivol,isid2,iobj2,count,refcnt,
     .        iobj_hold,isid_hold,inum,vwindex,iobj_primary,ipla,iint,
     .        problem_ignored,fp,ipro,iobj_hack,iobj_hack2
      LOGICAL cont,ref_debug
      TYPE(type_view) :: chord, tmpchord

      INTEGER, PARAMETER :: MAXREFNIND = 50000
c      INTEGER, PARAMETER :: MAXREFNIND = 5000000
      INTEGER reflvl,refnind,refind,reford,refend,refmark
      INTEGER, ALLOCATABLE :: refobj(:),refsid(:),refsrf(:),ref_hack(:)
      TYPE(type_view), ALLOCATABLE :: ref_chord(:)

      REAL   quant
      REAL*8 val,v1(3),v3(3),v4(3),v1_hold(3),v2_hold(3),s,length,frac

c      DATA problem_ignored /0/
      LOGICAL problem_message
      DATA    problem_message / .FALSE. /
      SAVE

      ref_debug = .FALSE.

      problem_ignored = 0

      output = .FALSE.
      
c     reflection chords are killed if this is activated
      dchord = -1 ! 78067 ! 8692 ! 7001 ! -1 ! -1 ! 6625  ! -1 
      fp = 0

      reford  = opt%ref_opt
      reflvl  = 1
      refind  = 0
      refmark = 0
      refnind = 0

c      IF (debug+1) THEN
c        WRITE(0,*) 'TERMINATING TRACE AFTER DIAGNOSTIC CHORD'
c        status = -2
c        RETURN
c      ENDIF

      chord = chord_primary
      chord%track => chord_primary%track
      chord%profile => chord_primary%profile
c      chord%spectrum => chord_primary%spectrum

      refcnt = 0
      iobj_hold = 1

c...  Decide if the detector is inside or outside the vessel wall:
      IF (opt%ccd.EQ.2) THEN
        vwindex = 1  ! C-Mod
      ELSE
        vwindex = 2
      ENDIF

c      WRITE(0,*) '------vwindex----',vwindex,opt%ccd  

      iobj_hack  = -1
      iobj_hack2 = -1

      cont = .TRUE.
      DO WHILE (cont)

        IF (nchord.EQ.dchord+1) THEN
          WRITE(0,*) 'TERMINATING TRACE AFTER DIAGNOSTIC CHORD'
          status = -2
          RETURN
        ENDIF

        cont = .FALSE.
c...    Find wall intersections for the chord, building VWINTER list of intersections:
        CALL FindSurfaceIntersections
     .         (chord%v1,chord%v2,IT_VWINTER,status)                       
c        WRITE(0,*) 'WALL INTESRSECTIONS:',nvwinter,nvwlist,vwindex
c        WRITE(0,'(10F12.6)') vwinter(1:nvwinter)%dist 
c        WRITE(0,*) vwinter(1:nvwinter)%obj 
c        WRITE(0,*) vwinter(1:nvwinter)%sur 
c        WRITE(0,*) obj(vwinter(1:nvwinter)%obj)%nsur
        IF (status.EQ.-1) THEN
          WRITE(0,'(A,10I8)') 
     .       'WARNING: NO WALL INTERSECTION FOUND'
          IF (refcnt.GT.0) THEN
            IF (iobj.GE.1.AND.iobj.LE.nobj) 
     .        WRITE(0,'(A,10I8)') '       DATA:',
     .          nchord,refcnt,iobj,obj(iobj)%ik,obj(iobj)%ir,count,
     .          iobj_primary,obj(iobj_primary)%ik,nvwlist
          ENDIF
          nchord = 1
          s_chord(nchord)%v1(1:3) = chord%v1(1:3)
          s_chord(nchord)%v2(1:3) = chord%v2(1:3)
          IF (refcnt.GT.0) THEN
            nchord = 2
            s_chord(nchord)%v1(1:3) = chord_primary%v1(1:3)
            s_chord(nchord)%v2(1:3) = chord_primary%v2(1:3)
          ENDIF
          status = -1
          RETURN
        ENDIF
c        WRITE(0,*) 'NCHORD, DCHORD:',nchord, dchord
        IF (output) THEN
          WRITE(fp,*) 'WALL INTERSECTIONS:',nvwinter,nvwlist
          WRITE(fp,*) vwinter(1:nvwinter)%dist 
          WRITE(fp,*) vwinter(1:nvwinter)%obj 
          WRITE(fp,*) vwinter(1:nvwinter)%sur 
          WRITE(fp,*) obj(vwinter(1:nvwinter)%obj)%nsur
          WRITE(fp,*) obj(vwinter(1:nvwinter)%obj)%index
        ENDIF
        IF ((refcnt.EQ.0.AND.nvwinter.LT.vwindex).OR.
     .      (refcnt.GT.0.AND.nvwinter.LT.1)) THEN
          WRITE(0,*) 'CHORD DOES NOT PASS CORRECTLY THROUGH VESSEL'
          WRITE(0,*) '  REFCNT=',refcnt,nvwinter,vwindex
          EXIT
        ELSEIF (refcnt.EQ.0) THEN
c...      Trim the end of the chord to where the line-of-sight
c         intersects the vessel wall:
          iobj = vwinter(vwindex)%obj
          iobj_primary = iobj
          chord%v2(1:3) = vwinter(vwindex)%v(1:3)  
        ELSEIF (refcnt.GT.0) THEN
          IF (vwinter(1)%dist.LT.1.0D-5) THEN
            refend = 2
          ELSE
            refend = 1
          ENDIF
          iobj          = vwinter(refend)%obj
          chord%v2(1:3) = vwinter(refend)%v(1:3) 
c          write(fp,*) 'reflected iobj...',
c     .      nchord,iobj,refend,vwinter(1)%dist
          tmpchord = chord
        ELSE
          CALL ER('IntegrateChord','Bad',*99)
        ENDIF

        v1_hold(1:3) = chord%v1(1:3)
        v2_hold(1:3) = chord%v2(1:3)

        IF (obj(iobj)%type.EQ.OP_INTEGRATION_VOLUME) THEN
c...      Extend the view a little to avoid precision problems when
c         this chord terminates inside the integration mesh:
c          WRITE(0,*) 'HERE!?',nchord
          length = DSQRT((chord%v2(1)-chord%v1(1))**2 + 
     .                   (chord%v2(2)-chord%v1(2))**2 + 
     .                   (chord%v2(3)-chord%v1(3))**2)
          frac = (length + 1.0D-7) / length
          chord%v2 = frac * (chord%v2 - chord%v1) + chord%v1


c   THETA=  1.664318555581006E-007
c   AV   =   7.069255      -2.422884      1.1765493E-06
c   V2   =   7.069731      -2.424536      1.1766284E-06
c   AV   = -4.7512795E-04  1.6522347E-03 -7.9076426E-11
c   BV   =  1.1766284E-06  0.0000000E+00  -7.069731    
c THETA:   129.376813239171     
c vv:  0.3960797      0.9182161      8.6442228E-08
c nv: -0.9610522     -0.2763668     -1.5994971E-07
c nv:  0.9610522      0.2763668      1.5994971E-07
c
c

        ENDIF

c...    Check for chord weight adjustment based on the type of surface being intersected:
        inum = obj(iobj)%index
        SELECTCASE (opt%obj_fudge(inum)) 
          CASE (0)
c           Do nothing:
          CASE (1)
            IF (refcnt.EQ.0) 
     .        chord%weight = chord%weight * opt%obj_factor(inum)
          CASEDEFAULT
            STOP 'UNRECOGNIZED FUDGE'
        ENDSELECT
  


c...    Find integration grid boundary intersections for the chord:
        IF (.TRUE.) THEN

c      WRITE(0,*) 'CHORD LENGTH: ',
c     .  DSQRT((chord%v2(1)-chord%v1(1))**2 + 
c     .        (chord%v2(2)-chord%v1(2))**2 + 
c     .        (chord%v2(3)-chord%v1(3))**2)

          ngbinter = 0
          CALL FindSurfaceIntersections
     .           (chord%v1,chord%v2,IT_GBINTER,status)                       

        ELSE
c...      Cheat, try last set of grid intersections first (if there is an appropriate one): 
        ENDIF

        IF (output) THEN
          WRITE(fp,*) 'GRID BOUNDARY INTERSECTIONS:',ngbinter
          WRITE(fp,*) gbinter(1:ngbinter)%dist 
          WRITE(fp,*) gbinter(1:ngbinter)%obj 
          WRITE(fp,*) gbinter(1:ngbinter)%sur 
          WRITE(fp,*) obj(gbinter(1:ngbinter)%obj)%nsur
        ENDIF

c        WRITE(fp,*) 'GRID BOUNDARY INTESRSECTIONS:',ngbinter
c        WRITE(fp,*) gbinter(1:ngbinter)%dist 
c        WRITE(fp,*) gbinter(1:ngbinter)%obj 
c        WRITE(fp,*) gbinter(1:ngbinter)%sur 
c        WRITE(fp,*) obj(gbinter(1:ngbinter)%obj)%nsur

c        IF (refcnt.GT.0) WRITE(0,*) 'GRID INTESRSECTIONS:',ngbinter


c...    So, this is necessary in case a chord is reflected from a wall surface
c       that's inside the integration volume, i.e. need to start the chord
c       off inside that volume rather than assuming it starts outside the grid: - 15/07/11, SL (was painful)
        iobj_hack2 = -1
        IF (refcnt.GT.0.AND.iobj_hack.NE.-1) THEN
c          WRITE(fp,*) 'ammending with hack....'
          IF (output) 
     .      WRITE(fp,*) 'ammending with hack....'
          IF (ref_debug) WRITE(6,*) 'ammending with hack....'
          gbinter(2:ngbinter+1) = gbinter(1:ngbinter)
          gbinter(1)%dist = 0.1D0 * gbinter(2)%dist  ! arbitrary
          gbinter(1)%obj  = iobj_hack
          gbinter(1)%sur  = 0
          gbinter(1)%v    = chord%v1(1:3)
          ngbinter = ngbinter + 1

c          WRITE(fp,*) 'GRID BOUNDARY INTESRSECTIONS:',ngbinter
c          WRITE(fp,*) gbinter(1:ngbinter)%dist 
c          WRITE(fp,*) gbinter(1:ngbinter)%obj 
c          WRITE(fp,*) gbinter(1:ngbinter)%sur 
c          WRITE(fp,*) obj(gbinter(1:ngbinter)%obj)%nsur

        ENDIF

        IF (ngbinter.GT.0) THEN  
c...      Intersections detected (NGBINTER.GT.0):
          DO i1 = 1, ngbinter, 2    ! The 2 since exit for each entrance, unless a floater is hit...

c        IF (output) THEN
c          WRITE(fp,*) 'check:',ngbinter,refcnt,
c     .       iobj_hold,obj(iobj_hold)%type.EQ.OP_INTEGRATION_VOLUME
c          WRITE(fp,*) '     :',SNGL(chord%v1(1:3))
c          WRITE(fp,*) '     :',iobj_hack
c        ENDIF

            IF (refcnt.GT.0.AND.
     .          obj(iobj_hold)%type.EQ.OP_INTEGRATION_VOLUME) THEN
c...          Check needed in case the reflection originates from inside
c             the integration mesh, i.e. from a target surface:
              iobj = iobj_hold
              isid = isid_hold
              v1(1:3) = chord%v1(1:3)
c            ELSEIF (refcnt.GT.0.AND.iobj_hack.NE.-1) THEN
c              iobj = iobj_hack
c              isid = 0
c              v1(1:3) = chord%v1(1:3)
            ELSE
              iobj = gbinter(i1)%obj
              isid = gbinter(i1)%sur
              v1(1:3) = gbinter(i1)%v(1:3)
            ENDIF

c            obj(iobj)%flag(isid) = -1         ! *** CHEAT ***

            count = 0
            DO WHILE (.TRUE.) 
              count = count + 1
c...          Assemble list of surfaces bounding this object: 
              noblist = 0
              DO i2 = 1, MAX(obj(iobj)%nsur,obj(iobj)%nside)     
c...            If the object geometry is toroidally discretized, then 
c               do not allow the surface that the chord is currently
c               "on" to be included in the search list 'cause ..., but otherwise okay...:
                IF (obj(iobj)%gsur(i2).EQ.GT_TD.AND.i2.EQ.isid) CYCLE     
c...            Add map surfaces, which were identified when objects were defined, to list: 
                IF (output) THEN
                  WRITE(fp,*) ' -->',iobj,i2,obj(iobj)%nmap(i2)
                ENDIF
                DO i3 = 1, obj(iobj)%nmap(i2) 
                  noblist = noblist + 1
                  oblist(noblist,1) = obj(iobj)%imap(i3,i2)   ! Need to store this list for the next ray... 
                  oblist(noblist,2) = obj(iobj)%isur(i3,i2)
                ENDDO
              ENDDO
              nobinter = 0

              IF (output) THEN
                WRITE(fp,*) ' ???:',iobj,isid
                WRITE(fp,*) '    :',obj(iobj)%ik,obj(iobj)%ir
                WRITE(fp,*) '    :',
     .            obj(iobj)%tsur(1:MAX(obj(iobj)%nsur,obj(iobj)%nside))
                WRITE(fp,*) '    :',noblist
                IF (noblist.GT.0) THEN
                  WRITE(fp,*) '    :',oblist(1:noblist,1)
                  WRITE(fp,*) '    :',obj(oblist(1:noblist,1))%ik
                  WRITE(fp,*) '    :',obj(oblist(1:noblist,1))%ir
                  WRITE(fp,*) '    :',oblist(1:noblist,2)
                ENDIF
              ENDIF


              CALL FindSurfaceIntersections
     .               (v1,chord%v2,IT_OBINTER,status)
c     .               (v1,chord%v2,IT_OBINTER,nobj,obj,status)

c              IF (refcnt.EQ.149) THEN
                IF (output) THEN
                  WRITE(fp,*) ' INT:',nobinter
                ENDIF
c              ENDIF



              IF (nobinter.EQ.0) THEN
c...            No intersections were found, check if a vessel wall
c               surface inside the integration volume was intersected:
            
c                WRITE(0,'(A,10I8)') 
c     .             'WARNING: NO OBJECT INTERSECTION FOUND',
c     .             nchord,refcnt,iobj,obj(iobj)%ik,obj(iobj)%ir,count
c                nchord = 1
c                s_chord(nchord)%v1(1:3) = v1_hold(1:3)
c                s_chord(nchord)%v2(1:3) = v2_hold(1:3)
c                status = -1
c                RETURN

c                v3(1:2) = chord%v2(1:2)            
c                v3(3)   = -100.0D0
c                v4(1:2) = v3(1:2)
c                v4(3)   = +100.0D0

c                ***THIS IS A WEAK CONDITION, BETTER TO CHECK IF END OF CHORD IS INSIDE
c                  CURRENT 3D VOLUME OF INTEREST, ALTHOUGH THIS MAY NOT BE STRAIGHT FORWARD
c                  SINCE THE VOLUME SOMETIMES HAS 'EXTRA' SURFACES ASSOCIATED WITH IT,
c                  FROM, LIKE FOR THE X-POINT WITH THE MAGNETIC GRID (NEEDS TO BE RESOLVED)***
                IF (.TRUE.) THEN

c...              Check if the end point of the chord is inside the current integration volume:
c                  IF (obj(iobj)%gsur(1).EQ.GT_TC) THEN
c                    WRITE(0,*) 'better',iobj
c                  ENDIF

c                  IF (i1.NE.ngbinter) problem_ignored =problem_ignored+1
                 problem_ignored = problem_ignored + 1
c                  IF (i1.NE.ngbinter) 
                  WRITE(0,*) 'PROBLEM IG...',nchord,refcnt,iobj
                  WRITE(0,*) '             ',obj(iobj)%ik,obj(iobj)%ir
                  WRITE(0,*) '             ',chord%v1(1:3)

                 IF (ref_debug) 
     .             write(6,*) 'refcnt,reflvl=',refcnt,reflvl
                 IF (refcnt.EQ.0) iobj_hack  = iobj
                 IF (refcnt.GT.0) iobj_hack2 = iobj
                 IF (ref_debug) 
     .             write(6,*) 'setting iobj_hack',iobj,iobj_hack,
     .                         iobj_hack2


c                IF (i1.EQ.ngbinter) THEN
c                IF (nobinter.EQ.2) THEN
c...              Put chord on vessel wall surface:
                  nobinter = 1
                  IF (refcnt.EQ.0) THEN 
                    obinter(1)%obj    = vwinter(vwindex)%obj
                    obinter(1)%sur    = vwinter(vwindex)%sur 
                    obinter(1)%v(1:3) = vwinter(vwindex)%v(1:3)   ! =chord%v2(1:3)
                  ELSE
                    obinter(1)%obj    = vwinter(refend)%obj
                    obinter(1)%sur    = vwinter(refend)%sur 
                    obinter(1)%v(1:3) = vwinter(refend)%v(1:3)   ! =chord%v2(1:3)
                  ENDIF
c                  obinter(1)%obj    = vwinter(1)%obj
c                  obinter(1)%sur    = vwinter(1)%sur 
c                  obinter(1)%v(1:3) = vwinter(1)%v(1:3)   ! =chord%v2(1:3)
                  obinter(1)%dist   = DSQRT((v1(1) - v2_hold(1))**2 +
     .                                      (v1(2) - v2_hold(2))**2 +
     .                                      (v1(3) - v2_hold(3))**2)

c                  obinter(1)%dist   = DSQRT((v1(1) - chord%v2(1))**2 +
c     .                                      (v1(2) - chord%v2(2))**2 +
c     .                                      (v1(3) - chord%v2(3))**2)
                  status = 0

                ELSE
c...              Problems:
                  STOP 'PROBLEMS...'
                ENDIF

              ELSEIF (nobinter.GE.2) THEN
c...            Try to identify when an intersection is invalid for a GT_TC surface:
                i2 = 1
                iobj2 = obinter(i2)%obj
                isid2 = obinter(i2)%sur
                IF (obj(iobj2)%gsur(isid2).EQ.GT_TC.AND.
     .              obinter(i2)%dist.LT.RTOL) THEN
c...              Check if the ISID2 surface is the same as the ISID (the surface
c                 that the chord is currently on), and if so, delete the intersection:
                  IF (isid.EQ.0) THEN
                    write(0,*) 'strange event here, ISID=0'
                    write(6,*) 'strange event here, ISID=0'
                    obinter(i2)%dist = -999.0D0
                  ELSE
                    DO i3 = 1, obj(iobj)%nmap(isid)
                      IF (obj(iobj)%imap(i3,isid).EQ.iobj2.AND.
     .                    obj(iobj)%isur(i3,isid).EQ.isid2)
     .                  obinter(i2)%dist = -999.0D0
                    ENDDO
                  ENDIF
                ENDIF
c...            Clean-up the list:
                IF (obinter(i2)%dist.EQ.-999.0D0) THEN
                  DO i3 = i2, nobinter-1
                    obinter(i3) = obinter(i3+1)
                  ENDDO
                  nobinter = nobinter - 1
                ENDIF
              ENDIF
 
              IF (nchord.EQ.-1) THEN
                WRITE(fp,'(A,5I6)') 
     .            ' DEBUG:',nchord,iobj,isid,obj(iobj)%ik,obj(iobj)%ir

                DO i2 = 1, nobinter
                  WRITE(fp,'(1X,A,F20.12,2I7)') '   OBINTER:',
     .              obinter(i2)%dist,
     .              obinter(i2)%obj,
     .              obinter(i2)%sur
                ENDDO
              ENDIF

              IF (status.LT.0) THEN
c...            Error condition:
                WRITE(0,*) '---PROBLEM WITH CHORD---'
                WRITE(0,*) '  STATUS   :',status
                WRITE(0,*) '  XYZ      :',chord%v1(1:3)
                WRITE(0,*) '  ANG      :',chord%xangle,chord%yangle
                WRITE(0,*) '  IOBJ,ISID:',iobj,isid
                WRITE(0,*) '  IK,IR    :',obj(iobj)%ik,obj(iobj)%ir
                DO i2 = 1, MAX(obj(iobj)%nsur,obj(iobj)%nside)
                  IF (i2.EQ.isid) CYCLE     
                  WRITE(0,'(A,I3,10(I7,2I4,2X))') 
     .              '   MAP:',i2,(obj(iobj)%imap(i3,i2),
     .                         obj(obj(iobj)%imap(i3,i2))%ik,
     .                         obj(obj(iobj)%imap(i3,i2))%ir,
     .                         i3=1,obj(iobj)%nmap(i2))
                  WRITE(0,'(A,I3,10(I7,8X,2X))') 
     .              '   SUR:',i2,(obj(iobj)%isur(i3,i2),
     .                         i3=1,obj(iobj)%nmap(i2))
                ENDDO
                WRITE(0,*) '  NGBINTER:',ngbinter
                WRITE(0,*) '  NOBINTER:',nobinter
                DO i2 = 1, nobinter
                  WRITE(0,'(1X,A,F10.6,2I7)') '          :',
     .              obinter(i2)%dist,
     .              obinter(i2)%obj,
     .              obinter(i2)%sur
                ENDDO               

                WRITE(0,*) '  HALTING INTEGRATION FOR CHORD',nchord

c...            *TEMP* (won't work with MPI...) 
                nchord = 1
                IF (nchord.LE.22500) THEN
                  s_chord(nchord)%v1(1:3) = v1(1:3)
                  s_chord(nchord)%v2(1:3) = chord%v2(1:3)
                ENDIF
                RETURN

              ELSEIF (count.GT.22500) THEN
                WRITE(0,*) '  COUNT LIMIT EXCEEDED (22500), HALTING '//
     .                     'INTEGRTATION FOR CHORD',nchord
                status = -1
                nchord = 1
                IF (nchord.LE.22500) THEN
                  s_chord(nchord)%v1(1:3) = v1(1:3)
                  s_chord(nchord)%v2(1:3) = chord%v2(1:3)
                ENDIF
                RETURN
              ENDIF
 
              IF (obj(iobj)%type.EQ.OP_INTEGRATION_VOLUME) THEN  ! Aren't all objects integration volumes here?  Or do floating surfaces pass as well?
c...            Line-of-sight integral:

                DO iint = 1, MAX(1,opt%int_num)
                  val = chord%weight * obinter(1)%dist * 
     .                  DBLE(obj(iobj)%quantity(iint))
                  chord%integral(iint) = chord%integral(iint) + val
               
c     .                                 chord%weight * 
c     .                                 obinter(1)%dist * 
c     .                                 DBLE(obj(iobj)%quantity(iint))

                IF (ref_debug)
     .             write(6,'(A,2I6,F10.4,3F10.4,1P,E10.2,0P,
     .                      I6,1P,E10.2,0P)') 
     .              '   weight',nchord,refcnt,chord%weight,
     .              obinter(1)%v(1:3),chord%integral(iint)/(4*3.14),
     .              iobj,obj(iobj)%quantity(iint)

                  IF (.FALSE..AND.refcnt.EQ.0) THEN   ! *** PROFILE HACK ***

c           WRITE(6,'(A,8F10.6)') ' distance:', 
c     .    chord_primary%v1(1:3),
c     .    obinter(1)%v(1:3),
c     .                DSQRT((obinter(1)%v(1) - chord_primary%v1(1))**2 +
c     .                      (obinter(1)%v(2) - chord_primary%v1(2))**2 +
c     .                      (obinter(1)%v(3) - chord_primary%v1(3))**2),
c     .            obinter(1)%dist

                    IF (iint.EQ.1) THEN
c                    IF (iint.EQ.1) THEN !.AND.obj(iobj)%quantity(1).GT.0.0) THEN
                      chord%nprofile = chord%nprofile + 1

                      ipro = chord%nprofile

                      chord%profile(ipro,-5) = 
     .                  DSQRT((obinter(1)%v(1)-chord_primary%v1(1))**2 +
     .                        (obinter(1)%v(2)-chord_primary%v1(2))**2 +
     .                        (obinter(1)%v(3)-chord_primary%v1(3))**2)-
     .                  0.5D0 * obinter(1)%dist
                      chord%profile(ipro,-4) = obinter(1)%dist
                      chord%profile(ipro,-3) = chord%weight
                      ipla = obj(iobj)%index_pla
                      chord%profile(ipro,-2) = plasma(ipla)%nb
                      chord%profile(ipro,-1) = plasma(ipla)%te
                      chord%profile(ipro, 0) = plasma(ipla)%tb

                      chord%profile(ipro,-7) = plasma(ipla)%nD
                      chord%profile(ipro,-6) = plasma(ipla)%nD2

                      chord%profile(ipro,-12) = DBLE(REAL(obj(iobj)%in))
                    ENDIF

                    IF (opt%int_charge(iint).GT.3) 
     .                STOP 'NOT ENOUGH PROFILE SPACE'
                    IF (opt%int_z(iint).GT.1)
     .                chord%profile(ipro,-11+opt%int_charge(iint)) =  
     .                   plasma(ipla)%ni(iint)                        

c                    write(0,*) 'test',-11+opt%int_charge(iint),
c     .                         plasma(ipla)%ni(iint)

                    chord%profile(ipro,iint) = 
     .                  chord%weight * DBLE(obj(iobj)%quantity(iint))
                  ENDIF


c...              Line shape:
                  IF (opt%int_type(iint).EQ.2) THEN
                    CALL CalculateLineShape                         ! I don't like this chord primary business!
     .                     (iint,chord,iobj,val,
     .                      obinter(1)%dist,v1,obinter(1)%v)
                  ENDIF
c...              Averaging:
                  IF (opt%int_type(iint).EQ.3) THEN
                    ipla = obj(iobj)%index_pla
                    SELECTCASE (opt%int_average(iint))
                      CASE (1) ! ne
                        quant = plasma(ipla)%ne
                      CASE (2) ! Te
                        quant = plasma(ipla)%te
                      CASE (3) ! nb
                        quant = plasma(ipla)%nb
                      CASE (4) ! vb 
                        quant = plasma(ipla)%vb
                      CASE (5) ! Tb
                        quant = plasma(ipla)%tb
                      CASE (6) ! ni_imp
                        quant = plasma(ipla)%ni(iint)
                      CASE (7) ! vi_imp
                        quant = plasma(ipla)%vi(iint)
                      CASE (8) ! Ti_imp
                        quant = plasma(ipla)%ti(iint)
                      CASE DEFAULT
                    ENDSELECT
                    chord%average(iint) = chord%average(iint) + 
     .                                    val * DBLE(quant)
                  ENDIF
                ENDDO
c...            Update inversion map based on track length in object volume: 
                ivol = obj(iobj)%ivolume
                chord%track(ivol) = chord%track(ivol) + 
     .                              chord%weight * 
     .                              obinter(1)%dist
c...            Keep track of sampling weight for object:
c                obj(ivol)%sample = obj(ivol)%sample + obinter(1)%dist  ! Done with %path... 
                IF (nchord.EQ.-1) WRITE(0,*) '  DIST:',obinter(1)%dist
              ENDIF
c              write(fp,*) '    iobj',iobj,chord%integral(1)  ! low check

c...          Update surface where start of view chord currently resides:
              iobj    = obinter(1)%obj
              isid    = obinter(1)%sur
              isrf    = obinter(1)%srf
              v1(1:3) = obinter(1)%v(1:3)

c              obj(iobj)%flag(isid) = -1  ! *** CHEAT ***

c...          Loop exit conditions:
              IF     (obj(iobj)%tsur(isid).EQ.SP_GRID_SURFACE) THEN
c...            Keep going:
              ELSEIF (obj(iobj)%tsur(isid).EQ.SP_GRID_BOUNDARY.OR.
     .                obj(iobj)%tsur(isid).EQ.SP_VESSEL_WALL) THEN
c...            Leave the loop:
                EXIT   
              ELSE
                WRITE(0,*) 'EXITING: *** PROBLEM WITH OBJECT ***',
     .            iobj,isid,refcnt,nchord
                WRITE(0,*) '                                    ',
     .            obj(iobj)%ik,obj(iobj)%ir,obj(iobj)%tsur(isid)
                status = -1
                nchord = 1
                IF (nchord.LE.22500) THEN
c                  s_chord(nchord)%v1(1:3) = tmpchord%v1(1:3)
c                  s_chord(nchord)%v2(1:3) = tmpchord%v2(1:3)
c                  nchord = nchord + 1
                  s_chord(nchord)%v1(1:3) = v1(1:3)
                  s_chord(nchord)%v2(1:3) = chord%v2(1:3)
                ENDIF
                RETURN
c                EXIT
              ENDIF

            ENDDO  ! Integration volume loop

          ENDDO  ! Integration grid boundary loop

        ELSE
c...      Line-of-sight for the viewing chord does not pass through the
c         integration volume:
        ENDIF

        IF (nchord.LE.22500) THEN
c        IF (nchord.LE.22500.AND.refcnt.EQ.0) THEN  ! reflection debug
          s_chord(nchord)%v1(1:3) = chord%v1(1:3)
          s_chord(nchord)%v2(1:3) = chord%v2(1:3)
c          WRITE(0,'(A,3F10.3)') 'ref3 chord%v1,2=',chord%v1
c          WRITE(0,'(A,3F10.3)') '               =',chord%v2
c          WRITE(0,*) 'nchord >>>',nchord
        ENDIF

        IF (n_pchord.LE.22500.AND.refcnt.EQ.0) THEN 
          p_chord(n_pchord)%v1(1:3) = chord%v1(1:3)
          p_chord(n_pchord)%v2(1:3) = chord%v2(1:3)
        ENDIF

c...    Spawn reflection chords, if necessary, then loop and add to integral for chord:
        IF (opt%ref_num.GT.0) THEN

          IF (ref_debug)
     .      write(6,*) 'spawning reflection chord'

50        refcnt = refcnt + 1
          IF (refcnt.EQ.1) THEN
            IF (reford.GT.1.AND.reflvl.GT.1) THEN 
c            IF (reford.EQ.2.AND.reflvl.EQ.2) THEN ! old
              refind = refind + 1
              iobj = refobj(refind)
              isid = refsid(refind)
              isrf = refsrf(refind)
              iobj_hack = ref_hack(refind)
c              iobj_hack = ref_hack(refnind)  ! bug! 19/09/2011
              chord%v1     = ref_chord(refind)%v1
              chord%v2     = ref_chord(refind)%v2
              chord%weight = ref_chord(refind)%weight_save
              IF (ref_debug) THEN
                WRITE(6,*) 'Pulling up a ghost...',reflvl,refind,
     .                     chord%weight
                WRITE(6,*) iobj,isid,isrf
                WRITE(6,*) iobj_hack
                WRITE(6,'(3F12.5)') chord%v1
                WRITE(6,'(3F12.5)') chord%v2
              ENDIF
            ELSE
              iobj = vwinter(vwindex)%obj
              isid = vwinter(vwindex)%sur
              isrf = vwinter(vwindex)%srf
            ENDIF
c...        Check if reflection model specified for this object side:
            if (ref_debug) write(6,*) obj(iobj)%reflec(isid)
            IF (obj(iobj)%reflec(isid).EQ.0) THEN
              refcnt = 0
            ELSE
              iobj_hold = iobj
              isid_hold = isid
            ENDIF
          ELSE
            iobj = -1
            isid = -1
            isrf = -1
c...        LEVEL=2: LEVEL=1 store the current reflection chord IOBJ,ISID,ISRF values for 
c           application later, as above, as well as the wall intersections, which are used in 
c           AssignReflectionChord. Set REFCNT to 1.
            IF (reford.GT.1.AND.reford.GT.reflvl) THEN
c            IF (reford.EQ.2.AND.reflvl.EQ.1) THEN  ! old
              IF (.NOT.ALLOCATED(ref_chord)) THEN
                ALLOCATE(refobj   (MAXREFNIND))                     !  wrong way to do things
                ALLOCATE(refsid   (MAXREFNIND))                     !  should just follow each chord to it's bitter end rather than
                ALLOCATE(refsrf   (MAXREFNIND))                     !  storing all the primary reflections in a list
                ALLOCATE(ref_hack (MAXREFNIND))
                ALLOCATE(ref_chord(MAXREFNIND))
              ENDIF
              refnind = refnind + 1
              IF (refnind.GT.MAXREFNIND) 
     .          CALL ER('IntegrateChord','Increase MAREFNIND',*99)
              refobj(refnind) = vwinter(refend)%obj
              refsid(refnind) = vwinter(refend)%sur
              refsrf(refnind) = vwinter(refend)%srf
              ref_hack(refnind) = iobj_hack2
              ref_chord(refnind) = chord
              IF (ref_debug) THEN
                WRITE(6,*) 'Storing reflected chord...',refnind,refend,
     .                      reflvl,chord%weight
                WRITE(6,*) 'iobj_hack=', iobj_hack
              ENDIF
            ENDIF
          ENDIF

          IF (refcnt.GT.0) 
     .       CALL AssignReflectionChord(chord,iobj,isid,isrf,
     .                                  refcnt,reflvl)          

c...      LEVEL=2:  Ugly, but if Assign... tries to trigger a loop exit, need to prep the 
c         next LEVEL=1 reflection chord to be processed need to spagetti back up.  Need to set
c         chord%weight to chord%weight*chord%weight

          IF (refcnt.EQ.0) THEN
            IF (reford.GT.1.AND.refind.LT.refnind) THEN
c            IF (reford.EQ.2.AND.refind.LT.refnind) THEN  ! old
              IF (ref_debug) 
     .          WRITE(6,*) 'Trying to leave...',refind,refnind

              IF (refind.EQ.refmark) THEN 
                reflvl  = reflvl + 1
                refmark = refnind
                IF (ref_debug) 
     .            WRITE(6,*) 'Moving up...',reflvl,refmark,refind
              ENDIF

              refcnt = 0
              GOTO 50
            ENDIF
          ELSE
            cont = .TRUE.
          ENDIF
c          IF (refcnt.NE.0) cont = .TRUE.

          IF (cont) nchord = nchord + 1  ! Debugging


c          IF (cont) THEN
c            tmpchord = chord
c            WRITE(0,*) 'REF:',chord%v1(1:3),refcnt
c            WRITE(0,*) '   :',chord%v2(1:3)
c          ENDIF

c reflection debug
c          IF (nchord.LE.22500.AND.refcnt.GT.0) THEN
c            s_chord(nchord)%v1(1:3) = chord%v1(1:3)
c            s_chord(nchord)%v2(1:3) = chord%v2(1:3)
c            WRITE(0,*) 'nchord,refcnt=',nchord,refcnt
c            WRITE(0,*) '   v1=',chord%v1(1:3)
c            WRITE(0,*) '   v2=',chord%v2(1:3)
c          ENDIF

        ENDIF

      ENDDO
      

c... *TEMP* (won't work with MPI...) 
c      IF (nchord.LE.22500) THEN
c        s_chord(nchord)%v1(1:3) = chord%v1(1:3)
c        s_chord(nchord)%v2(1:3) = chord%v2(1:3)
c      ENDIF
c      WRITE(0,*) 'SCHORD:',nchord
c      WRITE(0,*) 'SCHORD:',chord%v2(1:3)

c      WRITE(0,*) 'VIEW STORED'

c...  Kind of a silly way of doing this...? No, because there will be more than one chord for
c     each primary chord... still, ugly...  ** NO! JUST MODIFY v1,v2 FOR REFLECTIONS and WEIGHTS, AND 
c     GET RID OF _PRIMARY BUSINESS! 
      DO i1 = 1, MAX(1,opt%int_num)
        chord_primary%integral(i1) = chord%integral(i1)
        chord_primary%average (i1) = chord%average (i1)
      ENDDO
      chord_primary%nprofile = chord%nprofile

      write(0,*) 'iobj exit',iobj,isrf,obj(iobj)%iside(1,1:2),
     .             chord%weight
      srf(isrf)%count = srf(isrf)%count + chord%weight

      
c      DO i1 = -5, MAX(1,opt%int_num)  ! not required because pointers being used
c        chord_primary%profile(:,i1) = chord%profile(:,i1)
c        chord_primary%profile(:,2) = chord%profile(:,2)
c      ENDDO

      IF (problem_ignored.GE.1) THEN
        IF (.NOT.problem_message) THEN
          WRITE(0,*) 
          WRITE(0,*) '=============================================='
          WRITE(0,*) '"PROBLEM IGNORED"...',nchord,problem_ignored
          WRITE(0,*) '=============================================='
          WRITE(0,*) 
          problem_message = .TRUE.
        ENDIF
        IF (opt%ref_num.NE.0.AND.refcnt.GT.0) THEN
c        IF (.NOT.(opt%ref_num.NE.0.AND.refcnt.GT.0)) THEN
          STOP 'SHOULD NOT BE HERE RIGHT NOW'
          chord_primary%integral = 0.0
          chord_primary%average  = 0.0
        ENDIF
      ENDIF

c      WRITE(0,*) 'INTEGRAL STORED'

c      chord_primary%spectrum(1:100) = chord%spectrum(1:100)

c      WRITE(0,*) 'SPECTRUM STORED'

c      WRITE(0,*) 'DONE INTEGRATING CHORD+REFLEC'

      IF (ALLOCATED(ref_chord)) THEN
        DEALLOCATE(refobj)
        DEALLOCATE(refsid)
        DEALLOCATE(refsrf)
        DEALLOCATE(ref_hack)
        DEALLOCATE(ref_chord)
      ENDIF

c      WRITE(0,*) 'NCHORD, DCHORD:',nchord, dchord

      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE IntegrateAlongChords(pixel,status)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

c      TYPE(type_options985) :: opt
      TYPE(type_view)       :: pixel
c      TYPE(type_3D_object)  :: obj(nobj)
      INTEGER status                   
c      INTEGER nobj,status                   


      INTEGER i1,ix,iy,nxbin,nybin,n,iobj,iint,nview,save_ccd
      REAL*8  xangle,yangle,dxangle,dyangle,fact,mat(3,3),angle
c      TYPE(type_view) :: pixel2,chord
      TYPE(type_view) :: chord
 
c      INTEGER, ALLOCATABLE, TARGET :: itest(:)
c      INTEGER, POINTER :: ptr1(:)

      REAL*8, ALLOCATABLE, TARGET :: ddum1(:),ddum2(:,:)
      REAL*4, ALLOCATABLE, TARGET :: rdum1(:)


c      pixel2 = pixel_primary  ! Necessary?  Yes, so that the intergral value is passed back correctly?





c...  Object and pixel data is available.  Routine needs to be fast. 


c...  Want to process 1 pixel at a time on a particular processor (with
c     the option of chords for a given pixel originating at other pixels
c     if the focus is bad).  

c...  Build primary chords:

      pixel%integral = 0.0D0
      pixel%average  = 0.0D0

c      pixel2%integral = 0.0

      xangle  = pixel%xangle
      yangle  = pixel%yangle
      dxangle = pixel%dxangle
      dyangle = pixel%dyangle

      nxbin = pixel%nxbin
      nybin = pixel%nybin

c      xangle  = pixel2%xangle
c      yangle  = pixel2%yangle
c      dxangle = pixel2%dxangle
c      dyangle = pixel2%dyangle

c      nxbin = pixel2%nxbin
c      nybin = pixel2%nybin

      n = opt%n

      ALLOCATE(ddum1(n                          ))
      ALLOCATE(ddum2(nobj,-12:MAX(1,opt%int_num)))  ! MPI problem?  nobj=m, should be # integration volumes
      ALLOCATE(rdum1(100                        ))

      status = 0
      nview  = 0

      DO ix = 1, nxbin
        IF (status.LT.0) EXIT

        DO iy = 1, nybin
          IF (status.LT.0) EXIT

          status = 0

          nchord = nchord + 1
          n_pchord = n_pchord + 1
          nview  = nview  + 1

          ddum1 = 0.0D0
          ddum2 = 0.0D0

          chord%otrack = 1
          chord%track => ddum1
          chord%nprofile = 0
          chord%profile => ddum2  ! Waste of memory?  Just do a direct mapping from PIXEL...

          rdum1 = 0.0D0
c          chord%spectrum => rdum1

          chord%index_spe = pixel%index_spe
          chord%tmp       = pixel%tmp

          chord%weight = pixel%weight / DBLE(nxbin*nybin)
          chord%integral = 0.0D0
          chord%average  = 0.0D0
          chord%parent = pixel%index  ! Needed?

          DO i1 = 1, 3
c...        Perhaps add a small shift eventually? ...no, shouldn't bother 
            chord%v1(i1)    = pixel%v1(i1)
            chord%rot(i1)   = pixel%rot(i1)
            chord%trans(i1) = pixel%trans(i1)
          ENDDO

          fact = 0.5D0 * (1.0D0 / DBLE(nxbin) - 1.0D0) +
     .           DBLE(ix - 1) / DBLE(nxbin)
          
          chord%xangle  = xangle + fact * dxangle 
          chord%dxangle = 0.0D0

          fact = 0.5D0 * (1.0D0 / DBLE(nybin) - 1.0D0) +
     .           DBLE(iy - 1) / DBLE(nybin)
          
          chord%yangle  = yangle + fact * dyangle 
          chord%dyangle = 0.0D0

c          WRITE(0,*) 'CHORD:',nchord,chord%xangle,chord%yangle     

          IF (pixel%type.EQ.2) THEN
            save_ccd = opt%ccd
            opt%ccd  = pixel%ccd
            chord%v2 = pixel%v2
            WRITE(0,*) '---- chord%v1 ---->',SNGL(chord%v1)
            WRITE(0,*) '---- chord%v2 ---->',SNGL(chord%v2)
          ELSE
            chord%v2(1) = DSIN(chord%xangle * D_DEGRAD) * 
     .                    DCOS(chord%yangle * D_DEGRAD) * 
     .                    50.0D0 + chord%v1(1)
            chord%v2(2) = DCOS(chord%xangle * D_DEGRAD) * 
     .                    DSIN(chord%yangle * D_DEGRAD) * 
     .                    50.0D0 + chord%v1(2)
            chord%v2(3) = DCOS(chord%xangle * D_DEGRAD) * 
     .                    DCOS(chord%yangle * D_DEGRAD) * 
     .                    50.0D0 + chord%v1(3)
c            chord%v2(3) = -chord%v2(3)
	   
c...        Rotate and translate chords according to the detector specifications:

c...        Rotate about z-axis (roll):
            CALL Calc_Transform2(mat,0.0D0,1,0)
            angle = chord%rot(3)
            CALL Calc_Transform2(mat,angle,3,1)
            CALL Transform_Vect(mat,chord%v1)
            CALL Transform_Vect(mat,chord%v2)
c...        Rotate about x-axis (tilt):                       !... better to do swing before tilt? 
            CALL Calc_Transform2(mat,0.0D0,1,0)
            angle = chord%rot(1)
            CALL Calc_Transform2(mat,angle,1,1)
            CALL Transform_Vect(mat,chord%v1)
            CALL Transform_Vect(mat,chord%v2)
c...        Rotate about y-axis (swing):
            CALL Calc_Transform2(mat,0.0D0,1,0)
c...        BUG? 
            angle = chord%rot(2)
            CALL Calc_Transform2(mat,angle,2,1)
            CALL Transform_Vect(mat,chord%v1)
            CALL Transform_Vect(mat,chord%v2)
c...        Translate:
            chord%v1(1:3) = chord%v1(1:3) + chord%trans(1:3)
            chord%v2(1:3) = chord%v2(1:3) + chord%trans(1:3)
c            chord%v1(1:2) = chord%v1(1:2) + chord%trans(1:2)  
c            chord%v2(1:2) = chord%v2(1:2) + chord%trans(1:2)
c...        PHI rotoation:        
            angle = chord%trans(3) * 3.141596D0 / 180.0D0  ! PHI rotation added 29/03/2010 -SL
            CALL Calc_Transform2(mat,0.0D0,1,0)            ! Rotate about y-axis
            CALL Calc_Transform2(mat,angle,2,1)
            CALL Transform_Vect(mat,chord%v1)
            CALL Transform_Vect(mat,chord%v2)
          ENDIF

          IF (.FALSE.) THEN
c...        Toroidal rotation of camera, while preserving the view orientation with 
c           respect to the centre of the torus:
c            angle = 14.76D0 * 3.141596D0 / 180.0D0 ! LWIR
            angle = 18.7D0 * 3.141596D0 / 180.0D0 ! LWIR / MWIR?
c            angle = 4.04D0 * 3.141596D0 / 180.0D0
c            angle = 9.0D0 * 3.141596D0 / 180.0D0  ! For viewing one have of the plasma from HU12...
c...        Rotate about y-axis (swing):
            CALL Calc_Transform2(mat,0.0D0,1,0)
            CALL Calc_Transform2(mat,angle,2,1)
            CALL Transform_Vect(mat,chord%v1)
            CALL Transform_Vect(mat,chord%v2)
          ENDIF
c... 
          CALL IntegrateChord(chord,status)   ! new name

          IF (pixel%type.EQ.2) THEN
            opt%ccd  = save_ccd
          ENDIF

          IF (status.GE.0) THEN
c          IF (status.NE.-1) THEN
c...        Add integral result to pixel value:
            DO i1 = 1, MAX(1,opt%int_num)
              pixel%integral(i1) = pixel%integral(i1)+chord%integral(i1)
              pixel%average (i1) = pixel%average (i1)+chord%average (i1)
            ENDDO
c...
            pixel%track(1:n) = pixel%track(1:n) + chord%track(1:n)
            pixel%nprofile = MAX(pixel%nprofile,chord%nprofile)
            DO i1 = -12, MAX(1,opt%int_num)
             pixel%profile(:,i1)=pixel%profile(:,i1)+chord%profile(:,i1)
            ENDDO

c          write(6,*) '  check 1',pixel%profile(1:pixel%nprofile,1)

c            pixel%spectrum(1:100) = pixel%spectrum(1:100) + 
c     .                              chord%spectrum(1:100)

          ELSE
            WRITE(0,*) '  STATUS.LT.0 FOR CHORD'
            WRITE(0,*) '    PIXEL:',-1
            WRITE(0,*) '    CHORD:',ix,iy
          ENDIF

c...      Store representative vector for this pixel:
c          WRITE(0,*) ix,MAX(1,nxbin/2+1),iy,MAX(1,nybin/2+1),nchord
c          IF (ix.EQ.MAX(1,nxbin/2+1).AND.iy.EQ.MAX(1,nybin/2+1).AND.
c     .        nchord.LT.SIZE(s_chord,1)) THEN
c            pixel%global_v1 = SNGL(s_chord(nchord)%v1)
c            pixel%global_v2 = SNGL(s_chord(nchord)%v2)
          IF (ix.EQ.MAX(1,nxbin/2+1).AND.iy.EQ.MAX(1,nybin/2+1).AND.
     .        n_pchord.LT.SIZE(p_chord,1)) THEN
            pixel%global_v1 = SNGL(p_chord(n_pchord)%v1)
            pixel%global_v2 = SNGL(p_chord(n_pchord)%v2)
c
c
c            pixel%global_v1 = SNGL(chord%v1)
c            pixel%global_v2 = SNGL(chord%v2)
          ENDIF

        ENDDO
      ENDDO


c...  Finish up weighted average:
      DO iint = 1, opt%int_num
        IF (opt%int_type(iint).NE.3) CYCLE
        IF (DABS(pixel%integral(iint)).GT.1.0D-10) 
     .    pixel%average(iint) = pixel%average (iint) / 
     .                          pixel%integral(iint)
        WRITE(0,*) 'AVERAGE2:',iint,pixel%average(iint)
      ENDDO

c...  This averaging is bogus of course since the profiles won't overlap,
c     but is fine if there's a single chord, or if the view isn't too wide:  ! *** TRUE? ***
      DO i1 = -11, MAX(1,opt%int_num)
        pixel%profile(:,i1) = pixel%profile(:,i1) / DBLE(nview)
      ENDDO
c      pixel%profile(:,2) = pixel%profile(:,2) / DBLE(nview)

c      write(6,*) '  nview  ',nview
c      write(6,*) '  check 2',pixel%profile(1:pixel%nprofile,1)

c...  Clear memory:
      IF (ALLOCATED(ddum1)) DEALLOCATE(ddum1)
      IF (ALLOCATED(ddum2)) DEALLOCATE(ddum2)
      IF (ALLOCATED(rdum1)) DEALLOCATE(rdum1)

      RETURN
 99   STOP
      END
