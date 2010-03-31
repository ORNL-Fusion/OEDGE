c     -*-Fortran-*-
c
c ====================================================================== 
c
      SUBROUTINE AssignEmissionData(MAXPIXEL,npixel,pixel)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

      INTEGER MAXPIXEL,npixel
      TYPE(type_view) :: pixel(MAXPIXEL)

      INTEGER GetNumberOfObjects

      INTEGER iint,ipla,i2,iobj,ipixel,max_ik,max_ir,max_in,count,ndat,
     .        iobject,i1,i3,p(4)
      REAL    wlngth,delta
      CHARACTER fname*512
      REAL, ALLOCATABLE :: osm(:,:),tdata(:)



c...  Check for the presence of standard (OSM) and triangle (EIRENE) grid elements in 3D model:
      max_ik = 0
      max_ir = 0
      max_in = 0
      DO iobj = 1, nobj
        SELECTCASE (obj(iobj)%subtype)
          CASE (OP_FLUID_GRID)        
             max_ik = MAX(max_ik,obj(iobj)%ik)
             max_ir = MAX(max_ir,obj(iobj)%ir)
          CASE (OP_EIRENE_GRID)
             max_ik = MAX(max_ik,obj(iobj)%ik)
             max_ir = MAX(max_ir,obj(iobj)%ir)
             max_in = MAX(max_in,obj(iobj)%in)
          CASE DEFAULT
        ENDSELECT
      ENDDO


c...  Allocate memory for spectrum and plasma data storage:
      nplasma = 0
      nspectrum = 0
      count = 0
      DO iint = 1, opt%int_num
        IF (opt%int_type(iint).EQ.2) count = count + 1
      ENDDO
      IF (count.GT.0) THEN
        DO ipixel = 1, npixel
          IF (pixel(ipixel)%valid) nspectrum = nspectrum + 1
        ENDDO
        DO iobj = 1, nobj
          IF (obj(iobj)%type.EQ.OP_INTEGRATION_VOLUME) 
     .      nplasma = nplasma + 1
        ENDDO
        ALLOCATE(spectrum(MAXSPECBIN,count*nspectrum))
        ALLOCATE(plasma  (           count*nplasma  ))
        spectrum = 0.0
c        WRITE(0,*) 'CON<NPLA:',count,nplasma,count*nplasma
      ENDIF

      nplasma = 0
      nspectrum = 0
      DO iint = 1, opt%int_num

c...    Assign/calculate quantity of interest:

        SELECTCASE (opt%int_database(iint))
c         --------------------------------------------------------------
          CASE(1:3)
c          IF (max_ik.GT.0.AND.max_ir.GT.0) THEN
c...        Fluid grid quantity:
            IF (.NOT.ALLOCATED(osm)) ALLOCATE(osm(max_ik,max_ir))
            CALL GetFluidGridEmission(iint,max_ik,max_ir,osm,wlngth)
c         --------------------------------------------------------------
          CASE(4)
c           IF (.FALSE..AND.max_in.GT.0) THEN
c...        Take data from the EIRENE transfer file:

c...        Check which name 
            fname = 'default'
            i2 = 0
            i3 = 0
            DO iobject = 1, opt%obj_num
              IF (opt%obj_type(iobject).EQ.6) THEN
                DO i1 = 1, LEN_TRIM(opt%obj_fname(iobject))
                  IF (opt%obj_fname(iobject)(i1:i1).EQ.'.') THEN
                    IF (i2.EQ.0) i2 = i1
                    i3 = i1
                  ENDIF
c                  WRITE(0,*) opt%obj_fname(iobject)(i1:i1),i1,i2,i3
                ENDDO
                fname = 'eirene'//opt%obj_fname(iobject)(i2:i3)//
     .                  'transfer'
                WRITE(0,*) 'fname:',fname(1:LEN_TRIM(fname))
                EXIT
              ENDIF
            ENDDO

            ndat = GetNumberOfObjects(fname(1:LEN_TRIM(fname)))
            ALLOCATE(tdata(ndat))

            p = [0,0,0,0]
            IF (opt%int_line(iint).EQ.'B_ALPHA')    p = [6,1,7 ,1]  ! Dalpha
            IF (opt%int_line(iint).EQ.'B_GAMMA')    p = [6,2,6 ,1]  ! Dgamma
            IF (opt%int_line(iint).EQ.'C0_DENSITY') p = [2,2,1 ,1]  ! C atom density
            IF (opt%int_line(iint).EQ.'D+_DENSITY') p = [1,1,10,0]  ! D+ density
            WRITE(0,*) 'EIRENE DATA PARAMETERS=',p
            IF (p(1).EQ.0) CALL ER('AssignEmissionData','Bad '//
     .                             'EIRENE data tag',*99)
            CALL LoadTriangleData(p(1),p(2),p(3),p(4),tdata,fname)  
c            CALL LoadTriangleData(6,1,7 ,1,tdata,fname)  ! Dalpha
            WRITE(0,*) 'DONE'
c         --------------------------------------------------------------
          CASE DEFAULT
            CALL ER('ray_AssignEmssionData','Invalid database '//
     .              'specification',*99)
        ENDSELECT

        opt%int_wlngth(iint) = wlngth

c...    Assign data to objects:
        DO iobj = 1, nobj
          IF (obj(iobj)%type.NE.OP_INTEGRATION_VOLUME) CYCLE

          SELECTCASE (obj(iobj)%subtype)
            CASE (OP_FLUID_GRID)
c...          Fluid grid:
              obj(iobj)%quantity(iint) = osm(obj(iobj)%ik,obj(iobj)%ir) 
            CASE (OP_EIRENE_GRID)
c...          Eirene grid:
              IF (ALLOCATED(osm)) THEN
               obj(iobj)%quantity(iint) = osm(obj(iobj)%ik,obj(iobj)%ir) 
              ELSE
               obj(iobj)%quantity(iint) = tdata(obj(iobj)%in)
              ENDIF
            CASE (OP_INVERSION_GRID)
c...          Inversion mesh:

            CASE DEFAULT
              CALL ER('AssignEmissionData','Unknown integration '//
     .                'sub-type',*99)
          ENDSELECT

        ENDDO


        IF (opt%int_type(iint).EQ.2) THEN

c...      Setup spectrum binning:          
          delta = opt%int_width(iint) / REAL(MAXSPECBIN)
          opt%int_wlngthbin(1,iint) = -0.5*opt%int_width(iint)+0.5*delta
c     .      wlngth - 0.5 * opt%int_width(iint) + 0.5 * delta
          DO i2 = 2, MAXSPECBIN
            opt%int_wlngthbin(i2,iint)= opt%int_wlngthbin(i2-1,iint) + 
     .                                delta
          ENDDO
c          WRITE(0,*) opt%int_wlngthbin(1:MAXSPECBIN,iint)

c...      Map SPECTRUM array to PIXEL array:
          DO ipixel = 1, npixel
            IF (pixel(ipixel)%valid) THEN
              nspectrum = nspectrum + 1
              pixel(ipixel)%index_spe(iint) = nspectrum
            ENDIF
          ENDDO

        ENDIF
          

c...    Additional information for generating line shapes and weighted
c       averages (velocities, temperatures, B-field, ...):
        IF (opt%int_type(iint).EQ.2.OR.opt%int_type(iint).EQ.3) THEN

          DO iobj = 1, nobj
            IF (obj(iobj)%type.NE.OP_INTEGRATION_VOLUME) CYCLE
            IF (obj(iobj)%index_pla.EQ.0) THEN
              nplasma = nplasma + 1
              ipla = nplasma
              obj(iobj)%index_pla = ipla
            ELSE
              ipla = obj(iobj)%index_pla 
            ENDIF
            CALL AssignPlasmaQuantities(ipla,iint,iobj)
          ENDDO

        ENDIF


      ENDDO

   
      IF (ALLOCATED(osm  )) DEALLOCATE(osm  )
      IF (ALLOCATED(tdata)) DEALLOCATE(tdata)


      RETURN
 99   STOP
      END
c
c
c ======================================================================
c
      SUBROUTINE CalculateLineShape(iint,chord,iobj,val,dist,v1,v2)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

c...  Input:
      TYPE(type_view) chord
      INTEGER iobj,iint
      REAL*8 val,dist,v1(3),v2(3)

      INTEGER ipla,ispe,i1
      REAL    delta,mass,sigma,profile(MAXSPECBIN),protot,
     .        vel(3),los(3),len,dot,shift(MAXSPECBIN),
     .        shiftpro(MAXSPECBIN),bin(MAXSPECBIN),dopshift



      profile = 0.0

      ipla = obj(iobj)%index_pla

c...  Check if there are any emitters present:
c      WRITE(0,*) 'WHAT:',ipla,iint,plasma(ipla)%ni(iint)
      IF (plasma(ipla)%ni(iint).EQ.0.0) RETURN




c...  Static variables:



c     ...local binning to allow convolution...?  how is this done?



c...  Doppler broadening:     
      IF (.TRUE.) THEN
c...    Model:
        SELECTCASE (1)
          CASE (1)
c...        Simple normal distribution:
            delta = opt%int_width(iint) / REAL(MAXSPECBIN)
            mass = REAL(opt%int_a(iint)) * 1.67E-27   ! Store these....?
            sigma = opt%int_wlngth(iint) * 
     .              SQRT( (1.602E-19 * plasma(ipla)%ti(iint)) / 
     .                    (mass * 9.0E+16))
            DO i1 = 1, MAXSPECBIN
              profile(i1) = 1 / (sigma * 3.141592) *
     .                      EXP( -opt%int_wlngthbin(i1,iint)**2 / 
     .                            (2.0 * sigma**2) )
            ENDDO

c...        Normalization:
            protot = 0.0
            DO i1 = 1, MAXSPECBIN
              protot = protot + profile(i1) * delta
            ENDDO
            profile(1:MAXSPECBIN) = profile(1: MAXSPECBIN) / protot

            
          CASE DEFAULT
        ENDSELECT
      ENDIF

c      WRITE(0,*) mass,sigma,plasma(ipla)%ti(iint)
c      WRITE(0,*) protot
c      WRITE(0,*) profile


c...  Stark:
c...  Zeeman:


c...  Convolution?



c...  Doppler shift:
      IF (.TRUE.) THEN

c...    Get velocity data relative to LOS:
        SELECTCASE (obj(iobj)%type) 
          CASE (GT_TC) 

c...        Need to break the tranjectory up into descrete segments and 
c           evolve the plasma velocity vector across the trajectory:





          CASE (GT_TD) 

            bin(1:MAXSPECBIN) = opt%int_wlngthbin(1:MAXSPECBIN,iint)
            vel(1:3) = plasma(ipla)%bfield(1:3) * plasma(ipla)%vi(iint)
       

            len      = DSQRT((v2(1) - v1(1))**2 + 
     .                       (v2(2) - v1(2))**2 + 
     .                       (v2(3) - v1(3))**2)
            los(1:3) = SNGL((v2(1:3) - v1(1:3)) / len)

c            dot = plasma(ipla)%bfield(1) * los(1) + 
c     .            plasma(ipla)%bfield(2) * los(2) + 
c     .            plasma(ipla)%bfield(3) * los(3)
            dot = vel(1) * los(1) + vel(2) * los(2) + vel(3) * los(3)

c...        -ve DOT means flow towards viewer along line-of-sight 
            dopshift = opt%int_wlngth(iint) * 
     .                ((3.0E+8 + dot) / 3.0E+8 - 1.0)
            shift(1:MAXSPECBIN) = bin(1:MAXSPECBIN) + dopshift 

c            shift(1:MAXSPECBIN) = bin(1:MAXSPECBIN) *       ! Is this the right way to do it?
c     .                            (3.0E+8 + dot) / 3.0E+8 

            shiftpro = profile
            CALL Fitter(MAXSPECBIN,shift,shiftpro,
     .                  MAXSPECBIN,bin  ,profile ,'SPLINE')

c...        Check:
c            protot = 0.0
c            DO i1 = 1, MAXSPECBIN
c              protot = protot + profile(i1) * delta
c            ENDDO
c            WRITE(0,*) 'DAT:',protot,dopshift

          CASE DEFAULT
        ENDSELECT

c...    Model:
        SELECTCASE (1)
          CASE (1)
          CASE DEFAULT
        ENDSELECT

      ENDIF






c...  Scale (emissivity+dist+weight):

      ispe = chord%index_spe(iint)  

c      WRITE(0,*) '???',chord%index,ispe

      spectrum(1:MAXSPECBIN,ispe) = spectrum(1:MAXSPECBIN,ispe) + 
     .                              profile (1:MAXSPECBIN) * val
     .                         


c      WRITE(0,*) 'sdfsd',val,plasma(ipla)%ti(iint)

      RETURN
 99   STOP
      END



