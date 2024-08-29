c     -*-Fortran-*-
c
c ======================================================================
c
c
      SUBROUTINE DumpShinKajita(title9)
      USE mod_interface
      USE mod_out985
      USE mod_out985_variables
      use mod_params
      use mod_slout
      use mod_comgra
      use mod_cgeom
      use mod_walls_com
      use mod_dynam2
      use mod_dynam3
      use mod_pindata
      use mod_slcom
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'slout'
c     INCLUDE 'comgra'
c     INCLUDE 'cgeom'
c     INCLUDE 'walls_com'
c     INCLUDE 'dynam2'
c     INCLUDE 'dynam3'
c     INCLUDE 'pindata'
c     INCLUDE 'slcom'

      CHARACTER, INTENT(IN) :: title9*(*)

      INTEGER   ik,ir,fp,ike,ierr,iint,iobj,ik_last,ir_last,i
      REAL      fact(100),rdum
      CHARACTER dummy*1024,file*1024,tag*64
      
      CALL ZA09AS(dummy(1:8))
      dummy(9:10) = dummy(1:2)  ! Switch to EU format
      dummy(1:2 ) = dummy(4:5)
      dummy(4:5 ) = dummy(9:10)
      dummy(9:10) = '  '
      CALL ZA08AS(dummy(11:18))
      CALL CASENAME(dummy(21:),ierr)

c     E = h v; v = c / l; E = h c / l 
c      

      DO iint = 1, MAX(1,opt%int_num)
        write(0,*) 'wlngth',iint,opt%int_wlngth(iint)
           !  m2 kg / s  m / s    m
        fact(iint) = 6.63E-34 * 3.0E+8 / (opt%int_wlngth(iint) * 1.0E-9)
      ENDDO

      fp = 99
      OPEN (UNIT=fp,FILE='skf.impurity_lines',ACCESS='SEQUENTIAL',
     .      STATUS='REPLACE')
      WRITE(fp,'(A)') '# Data file for Kajita-san: line emission '//
     .                'in [W m-3]'
      WRITE(fp,'(A)') '#    (based on a SOLPS format from A. Kukushkin)'
      WRITE(fp,'(A)') '#'
      WRITE(fp,'(A)') '# Title       : '//TRIM(title9)
      WRITE(fp,'(A)') '# Case        : '//TRIM(dummy(21:))
      WRITE(fp,'(A)') '# Date & time : '//TRIM(dummy(1:18))
      WRITE(fp,'(A)') '# Version     : 1.0'
      WRITE(fp,'(A)') '#'
      WRITE(fp,'(A)') '# area - the area of the fluid code cell '//
     .                'in the poloidal plane'
      WRITE(fp,'(A)') '# vol  - the toroidal volume of the cell, '//
     .                'i.e. area*2*PI*x'
      WRITE(fp,'(A)') '#'


      WRITE(fp,'(A)') '# atomic number'
      WRITE(fp,'(10X,56X,25I14)') 
     .  (opt%int_z     (iint),iint=1,MAX(1,opt%int_num))

      WRITE(fp,'(A)') '# (approximate) atomic mass [amu]'
      WRITE(fp,'(10X,56X,25I14)') 
     .  (opt%int_a     (iint),iint=1,MAX(1,opt%int_num))

      WRITE(fp,'(A)') '# charge state'
      WRITE(fp,'(10X,56X,25I14)') 
     .  (opt%int_charge(iint),iint=1,MAX(1,opt%int_num))

      WRITE(fp,'(A)') '# wavelength [nm]'
      WRITE(fp,'(10X,56X,25F14.1)') 
     .  (opt%int_wlngth(iint),iint=1,MAX(1,opt%int_num))

      WRITE(fp,'(2A5,25A14)') '#  ix','iy','x','y','area','vol',
     .                       ('signal',iint=1,MAX(1,opt%int_num))
      WRITE(fp,'(2A5,25A14)') ' ',' ','[m]','[m]','[m-2]','[m-3]',
     .                       ('[W m-3]',iint=1,MAX(1,opt%int_num))

      ik_last = -1
      ir_last = -1

      DO iobj = 1, nobj
        IF (obj(iobj)%type   .NE.OP_INTEGRATION_VOLUME) CYCLE
        IF (obj(iobj)%gsur(1).NE.GT_TC                ) CYCLE

        ik = obj(iobj)%ik            
        ir = obj(iobj)%ir

        IF (ik.EQ.ik_last.AND.ir.EQ.ir_last) CYCLE

        ik_last = ik
        ir_last = ir         

        WRITE(fp,'(2I5,1P,25E14.4,0P)') 
     .    ik,ir,rs(ik,ir),zs(ik,ir),kareas(ik,ir),kvols(ik,ir),
     .    (obj(iobj)%quantity(iint)*fact(iint),
     .     iint=1,MAX(1,opt%int_num))

c        DO iint = 1, MAX(1,opt%int_num)
c        CALL inPutData(obj(iobj)%quantity(iint),TRIM(tag),
c     .                 'ph m-3 s-1')            
      ENDDO

 
      CLOSE(fp)


      IF (.TRUE.) THEN 
c
c       2D emission data: 
c       
        WRITE(file,'(A)') 'nc.divimp_imp_emission'
        CALL inOpenInterface(TRIM(file),NC_WRITE)
        CALL inPutData(opt%int_num,'N_SIGNAL','NA')
        DO i = 1, opt%int_num
          CALL inPutData(opt%int_z       (i),'ATOMIC_NUMBER','NA')
          CALL inPutData(opt%int_a       (i),'ATOMIC_MASS','NA')
          CALL inPutData(opt%int_charge  (i),'CHARGE','NA')
          CALL inPutData(opt%int_database(i),'DATABASE','NA')
          CALL inPutData(opt%int_wlngth  (i),'WAVELENGTH','nm')
        ENDDO

        DO iobj = 1, nobj
          IF (obj(iobj)%type.NE.OP_INTEGRATION_VOLUME) CYCLE

          ik = obj(iobj)%ik 
          ir = obj(iobj)%ir
          IF (ir.LT.irsep ) ik = ik - 1
          ir = ir - 1                    ! TUBE is set to the OSM fluid grid system, where
          IF (ir.GT.irwall) ir = ir - 2  ! the boundary rings are not present

          CALL inPutData(ik,'POS' ,'NA')
          CALL inPutData(ir,'TUBE','NA')
          DO i = 1, opt%int_num
            WRITE(tag,'(A,I0.2)') 'SIGNAL_',i
            CALL inPutData(obj(iobj)%quantity(i),tag,'ph m-3 s-1')
          ENDDO
        ENDDO
        CALL inCloseInterface

      ENDIF
 
      RETURN
99    STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE DumpLineData
      USE mod_interface
      USE mod_out985
      USE mod_out985_variables
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_dynam2
      use mod_dynam3
      use mod_outcom
      use mod_diagvel
      use mod_reiser_com
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'pindata'
c     INCLUDE 'dynam2'
c     INCLUDE 'dynam3'
c     INCLUDE 'outcom'
c     INCLUDE 'diagvel'
c     INCLUDE 'reiser_com'

      INTEGER   iint,iobj,ik,ir,ishift
      CHARACTER tag*7

      write(0,*) 'nobj',nobj

      CALL inOpenInterface('nc.line_dump',NC_WRITE)   ! TRIM(file) would not work, compiler bug...
      DO iint = 1, MAX(1,opt%int_num)
        CALL inPutData(opt%int_wlngth(iint),'WLNGTH','nm')
        WRITE(tag,'(A,I0.2)') 'LINE_',iint
        DO iobj = 1, nobj
c          write(0,*) 'test 1',obj(iobj)%type   .NE.OP_INTEGRATION_VOLUME
c          write(0,*) 'test 2',obj(iobj)%gsur(1).NE.GT_TC
          IF (obj(iobj)%type   .NE.OP_INTEGRATION_VOLUME) CYCLE
          IF (obj(iobj)%gsur(1).NE.GT_TC                ) CYCLE
          IF (iint.EQ.1) THEN
            ik = obj(iobj)%ik            
            ir = obj(iobj)%ir
            ishift = 1
            IF (ir.GT.irtrap) ishift = ishift + 2
            CALL inPutData(ik           ,'CELL','none')            
            CALL inPutData(ir-ishift    ,'TUBE','none')            
            CALL inPutData(kss(ik,ir)   ,'KSS','m')            
            CALL inPutData(kps(ik,ir)   ,'KPS','m')            
          ENDIF
          CALL inPutData(obj(iobj)%quantity(iint),TRIM(tag),
     .                   'ph m-3 s-1')            
        ENDDO
      ENDDO     
      CALL inCloseInterface

      RETURN
99    STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE AssignPlasmaQuantities(ipla,iint,iobj)
      USE mod_out985
      USE mod_out985_variables
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_dynam2
      use mod_dynam3
      use mod_outcom
      use mod_diagvel
      use mod_reiser_com
      IMPLICIT none

      INTEGER, INTENT(IN) :: ipla,iint,iobj

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'pindata'
c     INCLUDE 'dynam2'
c     INCLUDE 'dynam3'
c     INCLUDE 'outcom'
c     INCLUDE 'diagvel'
c     INCLUDE 'reiser_com'

c      INCLUDE 'slcom'

      INTEGER ik,ir,id,iz,zv,av,za
      REAL*8  Bx,By,Bz,beta,deltax,deltay,x(2),y(2),Btot,Btor,Bpol,
     .        mat(3,3),vector(3)

c      WRITE(0,*) 'NPLASM:',nplasma

      SELECTCASE (obj(iobj)%subtype)

        CASE (OP_FLUID_GRID)  ! *** PROFILE HACK *** (the OP_EIRENE_GRID)
c        CASE (OP_FLUID_GRID,OP_EIRENE_GRID)  
c...      Fluid grid:
          ik = obj(iobj)%ik
          ir = obj(iobj)%ir

          plasma(ipla)%nD = pinatom(ik,ir)
          plasma(ipla)%nD2= pinmol (ik,ir)
          plasma(ipla)%ne = knes (ik,ir)
          plasma(ipla)%te = ktebs(ik,ir)   ! Don't keep assigning this over and over...
          plasma(ipla)%nb = knbs (ik,ir)
          plasma(ipla)%vb = kvhs (ik,ir) / qtim
          plasma(ipla)%tb = ktibs(ik,ir)

c...      B-field components:         
          id = korpg(ik,ir)
          x(1) = DBLE(0.5 * (rvertp(1,id) + rvertp(2,id)))
          y(1) = DBLE(0.5 * (zvertp(1,id) + zvertp(2,id)))   
          x(2) = DBLE(0.5 * (rvertp(3,id) + rvertp(4,id)))
          y(2) = DBLE(0.5 * (zvertp(3,id) + zvertp(4,id)))   
          deltax = (x(2) - x(1))
          deltay = (y(2) - y(1))
          Btot = 1.0D0
          Bpol = DBLE(bratio(ik,ir))                              ! B_pol / B_tot
          IF (DABS(deltay).LT.1.0D-10) THEN
            beta = 0.0D0
          ELSE
            beta = deltax / deltay
          ENDIF 
          Btor = DSQRT(Btot**2 - Bpol**2)                                     
          Bz = Btor
          By = Bpol * DSQRT(1.0D0/(1.0D0+beta**2)) * DSIGN(1.0D0,deltay)
          Bx = beta * By
          plasma(ipla)%bfield(0) = 2.0  ! Tesla, just for laughs...
          plasma(ipla)%bfield(1) = Bx
          plasma(ipla)%bfield(2) = By
          plasma(ipla)%bfield(3) = Bz

c          IF (ir.EQ.irsep) THEN
c            WRITE(0,*) 'BFIELD:',plasma(ipla)%bfield(1:3),
c     .        obj(iobj)%phi*180.0/3.141592
c          ENDIF

c...      Rotate about y-axis:
          CALL Calc_Transform2(mat,0.0D0,1,0)
          CALL Calc_Transform2(mat,DBLE(obj(iobj)%phi),2,1)
          vector(1:3) = DBLE(plasma(ipla)%bfield(1:3))
          CALL Transform_Vect(mat,vector)
          plasma(ipla)%bfield(1:3) = SNGL(vector(1:3))

c          IF (ir.EQ.irsep) THEN
c            WRITE(0,*) 'BFIELD:',plasma(ipla)%bfield(1:3)
c          ENDIF

          za = opt%int_z(iint) * 1000 + opt%int_a(iint)

          SELECTCASE (za)

            CASE (01002)  ! Deuterium 
              SELECTCASE (opt%int_charge(iint))
                CASE (0)
                  plasma(ipla)%ni(iint) = knbs(ik,ir)   ! Temporary proxy data... 
                  plasma(ipla)%vi(iint) = kvhs(ik,ir) / qtim
                  plasma(ipla)%ti(iint) = ktebs(ik,ir)
                CASE DEFAULT
                  STOP 'SORRY, NO D BROADENING'
              ENDSELECT
    
            CASE DEFAULT
c...          Check if impurity data is requested, and whether or not it is 
c             available:
              zv = opt%int_z(iint)
              av = opt%int_a(iint) 
              iz = opt%int_charge(iint)

              IF     (iz.GE.0.AND.iz.LE.nizs) THEN
c                IF (.NOT.debugv)  ! *** PROFILE HACK ***
c     .            CALL ER('AssignPlasmaQuantities','Impurity velocity'//
c     .                    ' data not present',*99)
                plasma(ipla)%ni(iint) = sdlims(ik,ir,iz) * absfac
c                plasma(ipla)%vi(iint) = sdvs  (ik,ir,iz)  ! *** PROFILE HACK ***
                plasma(ipla)%ti(iint) = sdts  (ik,ir,iz)


c                IF (ir.EQ.irsep)
c     .                WRITE(0,*) sdvs(ik,ir,iz),velavg(ik,ir,iz)

c                IF (ir.EQ.irsep) THEN
c                  WRITE(0,*) 'VEL:',iz,ik,ir,sdvs(ik,ir,iz)
c                ENDIF
              ELSEIF (iz.EQ.0.AND..FALSE.) THEN
              ELSE
                CALL ER('AssignPlasmaQuantities','Impurity data '//
     .                  'not available',*99)
              ENDIF

          ENDSELECT

        CASE (OP_EIRENE_GRID)  ! *** PROFILE HACK ***
cc...      Eirene grid:
c          STOP 'NOT READY: OP_EIRENE_GRID'
        CASE (OP_INVERSION_GRID)
c...      Inversion mesh:
          STOP 'NOT READY: OP_INVERSION_GRID'
        CASE DEFAULT
          CALL ER('AssignPlasmaQuantities','Unknown integration '//
     .            'sub-type',*99)
      ENDSELECT


      RETURN
 99   WRITE(0,*) 'INTEGRATION CHANNEL: ',iint
      WRITE(0,*) '    Z =',zv,'  A=',av,'  CHARGE=',iz
      WRITE(0,*) '    ZA=',za
      STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE GetFluidGridEmission(iint,ik,ir,osm,wlngth2)
      USE mod_out985
      USE mod_out985_variables
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_dynam3
      use mod_pindata
      use mod_outcom
      IMPLICIT none

      INTEGER, INTENT(IN) :: iint,ik,ir
      LOGICAL display_warning 
      REAL    osm(ik,ir),wlngth2,scale

      DATA display_warning /.TRUE./

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'dynam3'
c     INCLUDE 'pindata'
c     INCLUDE 'outcom'


      INTEGER za,iz,ik1,ir1
      REAL    plrpad(MAXNKS,MAXNRS)

      za = opt%int_z(iint) * 1000 + opt%int_a(iint)

      SELECT CASE (opt%int_database(iint))
c       ----------------------------------------------------------------
        CASE (1)
c...      PIN: 
          IF (opt%int_line(iint).EQ.'UNITY') THEN
            osm = 1.0
            wlngth2 = 1.0
          ELSE
            SELECTCASE (za)
              CASE (01002) ! Deuterium

                IF     (opt%int_line(iint).EQ.'B_ALPHA') THEN
                  osm(1:ik,1:ir) = pinline(1:ik,1:ir,6,H_BALPHA)
c                 *** TEMP *** get rid of bad D2+ data at high temperatures
                  DO ir1 = 1, ir ! nrs
                    DO ik1 = 1, nks(ir1)
                      IF (ktebs(ik1,ir1).GT.1.0E+3.AND.
     .                  pinline(ik1,ir1,4,H_BALPHA).GT.
     .                  pinline(ik1,ir1,1,H_BALPHA)*10.0) THEN
c     .                  pinline(ik1,ir1,4,H_BALPHA).GT.1.0E+20) THEN
                        pinline(ik1,ir1,4,H_BALPHA) = 0.0
                        osm(ik1,ir1)=SUM(pinline(ik1,ir1,1:3,H_BALPHA))+ 
     .                                   pinline(ik1,ir1,5  ,H_BALPHA)
                        IF (display_warning) THEN
                          display_warning = .FALSE.
                          WRITE(0,*)
                          WRITE(0,*) '===================='
                          WRITE(0,*) '  CLIPPING D_ALPHA  '
                          WRITE(0,*) '===================='
                          WRITE(0,*)
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDDO
                  IF (.NOT.display_warning) THEN
                    DO ir1 = 1, ir ! nrs
                      DO ik1 = 1, nks(ir)
c                       write(0,*) 'debug',ik1,ir1,nrs
                       osm(ik1,ir1) = SUM(pinline(ik1,ir1,1:5,H_BALPHA))
                      ENDDO
                    ENDDO
                  ENDIF

                  wlngth2 = 656.3  ! Air

                ELSEIF (opt%int_line(iint).EQ.'B_GAMMA') THEN
                  osm(1:ik,1:ir) = pinline(1:ik,1:ir,6,H_BGAMMA)
                  wlngth2 = 434.0  ! Air

                ENDIF

              CASE DEFAULT
c...            ???:
                IF     (opt%int_line(iint).EQ.'B_ALPHA') THEN
                  osm(1:ik,1:ir) = pinline(1:ik,1:ir,6,H_BALPHA)
                  wlngth2 = 656.3  ! Air

                ELSEIF (opt%int_line(iint).EQ.'B_GAMMA') THEN
                  osm(1:ik,1:ir) = pinline(1:ik,1:ir,6,H_BGAMMA)
                  wlngth2 = 434.0  ! Air

                ENDIF

            ENDSELECT
          ENDIF      
c       ----------------------------------------------------------------   
        CASE (2)
c...      ADAS: 

          SELECTCASE (za)
            CASE (01002) ! Deuterium
              CALL LDADAS(opt%int_z(iint),iz,    ! requested matches what's available...
     .                    opt%int_adasid(iint),
     .                    opt%int_adasyr(iint),
     .                    opt%int_adasex(iint),
     .                    opt%int_isele (iint),
     .                    opt%int_iselr (iint),
     .                    opt%int_iselx (iint),
     .                    plrpad,wlngth2,ircode) 

                scale = 1.0
c              STOP 'D ADAS NOT READY'
            CASE DEFAULT
c...          Check for impurity data:
c             CALL LDADAS(1,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
c     >                   plrpad,wlngth2,IRCODE)

              iz = opt%int_charge(iint)

c              STOP 'PROBLEM WITH ADAS.F THAT I FIXED AND REMOVED'

              IF (iz.GE.0.AND.iz.LE.nizs) THEN     ! Add some checks to make sure the impurity 
                CALL LDADAS(opt%int_z(iint),iz,    ! requested matches what's available...
     .                      opt%int_adasid(iint),
     .                      opt%int_adasyr(iint),
     .                      opt%int_adasex(iint),
     .                      opt%int_isele (iint),
     .                      opt%int_iselr (iint),
     .                      opt%int_iselx (iint),
     .                      plrpad,wlngth2,ircode) 

                scale = absfac
              ELSE
                CALL ER('GetFluidGridEmission','Specified charge '//
     .                  'state invalid',*99)
              ENDIF
          ENDSELECT

          IF (ircode.NE.0) 
     .      CALL ER('GetFluidGridEmission','IRCODE.NE.0',*99)

          wlngth2 = wlngth2 / 10.0  ! A to nm... 

          osm(1:ik,1:ir) = plrpad(1:ik,1:ir) * scale
          WRITE(0,*) 'WLNGTH:',wlngth2
c       ----------------------------------------------------------------
        CASE (3) 
c...      DIVIMP pre-calculated quantities (should move DALHPA and DGAMMA here):
          IF     (TRIM(opt%int_line(iint)).EQ.'PRAD') THEN
            SELECTCASE (opt%int_charge(iint))
              CASE (-1)  ! Sum over all charge states 
                osm = 0.0
                DO iz = 0, MIN(cion,nizs)
c                  WRITE(0,*) 'IZ!',iz,absfac
                  osm(1:ik,1:ir) = osm(1:ik,1:ir) + powls(1:ik,1:ir,iz)
                ENDDO
                IF (absfac.GT.0.0) osm = osm * absfac
                DO iz = 0, 1
c                  WRITE(0,*) 'IZ!',iz
                  osm(1:ik,1:ir) = osm(1:ik,1:ir) + hpowls(1:ik,1:ir,iz)
                ENDDO

              CASE DEFAULT
                CALL ER('GetFluidGridEmission','Unknown PRAD opt',*99)
            ENDSELECT
          ELSE
            CALL ER('GetFluidGridEmission','Unknown DIVIMP tag',*99)
          ENDIF
c       ----------------------------------------------------------------
        CASE DEFAULT
          STOP 'UNRECOGNIZED DATABASE'
c       ----------------------------------------------------------------
      ENDSELECT

      RETURN
 99   WRITE(0,*) '    IZ=',iz,'  NIZS=',nizs
      STOP
      END
c
c ====================================================================== 
c
      LOGICAL FUNCTION CheckInversionCell(mode,x,y)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_slcom
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'slcom'

      INTEGER mode
      REAL*8  x,y

      INTEGER nlist,count,ik,ir,id,i1,i2
      LOGICAL init
      REAL*8  xlist(2,1000),ylist(2,1000),t1,t2,
     .        rvp(5,MAXNKS*MAXNRS),zvp(5,MAXNKS*MAXNRS),
     .        arxp,azxp,rfrac,zfrac

      DATA init /.TRUE./

      SAVE

      CheckInversionCell = .TRUE.

      IF (init) THEN
c...    List of boundary points:
        WRITE(0,*) 'BUILDING LIST!'
        init = .FALSE.
        nlist = 0
        xlist = 0.0D0
        ylist = 0.0D0
        rvp = DBLE(rvertp)
        zvp = DBLE(zvertp)

        IF (.TRUE.) THEN
c...      Stretch grid slightly:
          arxp = DBLE(rxp)
          azxp = DBLE(zxp)
          rfrac = 1.10D0
          zfrac = 1.05D0
          DO ir = 2, nrs
            IF (idring(ir).EQ.BOUNDARY) CYCLE
            DO ik = 1, nks(ir)
              id = korpg(ik,ir)
              DO i1 = 1, nvertp(id)
c...            Radial:
                IF     (rvp(i1,id).LT.arxp) THEN
c                  WRITE(0,*) 'ADJUSTING:',zvp(i1,id)
                  rvp(i1,id) =  arxp - rfrac * (arxp - rvp(i1,id))
c                  WRITE(0,*) '         :',zvp(i1,id)
                ENDIF
c...            Vertical:
                IF     (zvp(i1,id).GT. azxp) THEN
                  zvp(i1,id) =  azxp + zfrac * (zvp(i1,id) - azxp)
                ELSEIF (zvp(i1,id).LT.-azxp) THEN
                  zvp(i1,id) =  azxp + zfrac * (zvp(i1,id) - azxp)
                ENDIF
              ENDDO
            ENDDO
          ENDDO

        ENDIF

c...    Core:
        ir = 2
        DO ik = 1, nks(ir)-1
          id = korpg(ik,ir)
          nlist = nlist + 1
          xlist(1,nlist) = rvp(1,id)
          ylist(1,nlist) = zvp(1,id)
          xlist(2,nlist) = rvp(4,id)
          ylist(2,nlist) = zvp(4,id)
        ENDDO
c...    Everywhere else, radially:
         
        WRITE(0,*)
        WRITE(0,*) '  ================================== '
        WRITE(0,*) '  GRID CLIPPING HACK: RING 32 FORCED '
        WRITE(0,*) '  ================================== '
        WRITE(0,*)
        DO ir = 3, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
           STOP 'OLD HACK'
          IF (ir.EQ.36) CYCLE                      ! *** HACK ***
          DO ik = 1, nks(ir)-1
            IF     (irins(ik,ir).EQ.irwall.OR.
     .              irins(ik,ir).EQ.irtrap) THEN
              id = korpg(ik,ir)
              nlist = nlist + 1
              xlist(1,nlist) = rvp(1,id)
              ylist(1,nlist) = zvp(1,id)
              xlist(2,nlist) = rvp(4,id)
              ylist(2,nlist) = zvp(4,id)
            ELSEIF (irouts(ik,ir).EQ.irwall.OR.
     .              irouts(ik,ir).EQ.irtrap.OR.
     .              ir.EQ.32) THEN                 ! *** HACK ***
              id = korpg(ik,ir)
              nlist = nlist + 1
              xlist(1,nlist) = rvp(2,id)
              ylist(1,nlist) = zvp(2,id)
              xlist(2,nlist) = rvp(3,id)
              ylist(2,nlist) = zvp(3,id)
            ENDIF
          ENDDO
        ENDDO
c...    Targets:
        DO ir = irsep, nrs
          ik = 1
          id = korpg(ik,ir)
          nlist = nlist + 1
          xlist(1,nlist) = rvp(1,id)
          ylist(1,nlist) = zvp(1,id)
          xlist(2,nlist) = rvp(2,id)
          ylist(2,nlist) = zvp(2,id)
          ik = nks(ir)
          id = korpg(ik,ir)
          nlist = nlist + 1
          xlist(1,nlist) = rvp(3,id)
          ylist(1,nlist) = zvp(3,id)
          xlist(2,nlist) = rvp(4,id)
          ylist(2,nlist) = zvp(4,id)
        ENDDO
      ENDIF

c...  Check if cell is inside standard mesh:
      count = 0
      DO i1 = 1, nlist
        CALL CalcInter(xlist(1,i1),ylist(1,i1),xlist(2,i1),ylist(2,i1),
     .                 x,y,x+100.0D0,y,t1,t2)
        IF (t1.GE.0.0D0.AND.t1.LT.1.0D0.AND.
     .      t2.GE.0.0D0.AND.t2.LT.1.0D0) THEN
          count = count + 1
c          WRITE(0,*) '->',t1,t2,i1
c          WRITE(0,*) '  ',xlist(1,i1),ylist(1,i1)
c          WRITE(0,*) '  ',xlist(2,i1),ylist(2,i1)
        ENDIF
      ENDDO
      IF (MOD(count,2).EQ.0) CheckInversionCell = .FALSE.
c      WRITE(0,*) 'CELL CHECK:',x,y,count,CheckInversionCell


      RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE GetSepDist(mode,x1,x2,distmin)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_dynam2
      use mod_slcom
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'pindata'
c     include 'dynam2'
c     INCLUDE 'slcom'

      INTEGER mode
      REAL*8  x1,x2,distmin

      INTEGER CalcPoint

      INTEGER nir,ir1(10),i1,ik,ir,id   
      REAL*8  a1,a2,b1,b2,t,r,z,dist

      nir = 1
      ir1(1) = irsep
      IF (nrs.EQ.65) THEN   ! Special for PSI CDN grid 
        nir = nir + 1
        ir1(nir) = 38
      ENDIF

      distmin = DBLE(HI)

      DO i1 = 1, nir
        ir = ir1(i1)

        DO ik = 1, nks(ir)
 
          SELECTCASE (mode)
            CASE (1)
              id = korpg(ik,ir)
              a1 = DBLE(rvertp(1,id))
              a2 = DBLE(zvertp(1,id))
              b1 = DBLE(rvertp(4,id))
              b2 = DBLE(zvertp(4,id))

              IF (CalcPoint(a1,a2,b1,b2,x1,x2,t).EQ.2) THEN
                r    = a1 + t * (b1 - a1)
                z    = a2 + t * (b2 - a2)
                dist = DSQRT((r - x1)**2 + (z - x2)**2)
                IF (dist.LT.distmin) distmin = dist

c                WRITE(0,*) 'DATA:',x1,x2,ik,ir,dist,distmin
              ENDIF

            CASE(2)
              r = DBLE(rs(ik,ir))
              z = DBLE(zs(ik,ir))
              dist = DSQRT((r - x1)**2 + (z - x2)**2)
              IF (dist.LT.distmin) distmin = dist

            CASE DEFAULT
              CALL ER('GetSepDist','Bad mode',*99)

          ENDSELECT

        ENDDO

      ENDDO



      RETURN
 99   STOP
      END
c
c ======================================================================
c
       FUNCTION fq(xx, yy) RESULT(fn_val)
       implicit none
       REAL*8, INTENT(IN)  :: xx, yy
       REAL*8              :: fn_val

       fn_val = ((xx + 2.*yy)/3.) ** 2
       RETURN
       END FUNCTION fq

       FUNCTION fx(xx, yy) RESULT(fn_val)
       implicit none
       REAL*8, INTENT(IN)  :: xx, yy
       REAL*8              :: fn_val

       fn_val = 2. * (xx + 2.*yy) / 9.
       RETURN
       END FUNCTION fx

       FUNCTION fy(xx, yy) RESULT(fn_val)
       implicit none
       REAL*8, INTENT(IN)  :: xx, yy
       REAL*8              :: fn_val

       fn_val = 4. * (xx + 2.*yy) / 9.
       RETURN
       END FUNCTION fy


       SUBROUTINE qs2test2
!                           QS2TEST

!   THIS PROGRAM TESTS THE SCATTERED DATA INTERPOLATION PACKAGE QSHEP2D
! BY PRINTING THE MAXIMUM ERRORS ASSOCIATED WITH INTERPOLATED VALUES AND
! GRADIENTS ON A 10 BY 10 UNIFORM GRID IN THE UNIT SQUARE.  THE DATA SET
! CONSISTS OF 36 NODES WITH DATA VALUES TAKEN FROM A QUADRATIC FUNCTION FOR
! WHICH THE METHOD IS EXACT.  THE RATIO OF MAXIMUM INTERPOLATION ERROR
! RELATIVE TO THE MACHINE PRECISION IS ALSO PRINTED.  THIS SHOULD BE O(1).
! THE INTERPOLATED VALUES FROM QS2VAL AND QS2GRD ARE COMPARED FOR AGREEMENT.

       USE qshep2d
       IMPLICIT NONE
       INTEGER :: i, ier, j, k, lcell(3,3), lnext(36)
       REAL*8    :: a(5,36), dx, dy, eps, eq, eqx, eqy, f(36), 
     .            p(10), px, py, q, q1, 
     .            qx, qy, rmax, rq, rsq(36), x(36), xmin, y(36), yk, 
     .            ymin

       REAL*8 fx,fy,fq

! QSHEP2 PARAMETERS AND LOGICAL UNIT FOR OUTPUT

       INTEGER, PARAMETER  :: lout = 6, n = 36, nq = 13, nr = 3, nw = 19

! GENERATE A 6 BY 6 GRID OF NODES IN THE UNIT SQUARE WITH THE NATURAL ORDERING.

       k = 0
       DO  j = 1, 6
         yk = (6-j) / 5.
         DO  i = 1, 6
           k = k + 1
           x(k) = (i-1) / 5.
           y(k) = yk
         END DO
       END DO

! COMPUTE THE DATA VALUES.

       DO  k = 1, n
         f(k) = fq(x(k),y(k))
       END DO

! COMPUTE PARAMETERS DEFINING THE INTERPOLANT Q.

       CALL qshep2(n,x,y,f,nq,nw,nr,lcell,lnext,xmin,ymin,
     .             dx,dy,rmax,rsq,a,ier)
       IF (ier == 0) THEN

! GENERATE A 10 BY 10 UNIFORM GRID OF INTERPOLATION POINTS (P(I),P(J)) IN THE
!   UNIT SQUARE.  THE FOUR CORNERS COINCIDE WITH NODES.

         DO  i = 1, 10
           p(i) = (i-1) / 9.
         END DO

! COMPUTE THE MACHINE PRECISION EPS.

         eps = EPSILON(1.0)

! COMPUTE INTERPOLATION ERRORS AND TEST FOR AGREEMENT IN THE
!   Q VALUES RETURNED BY QS2VAL AND QS2GRD.

         EQ = 0.
         eqx = 0.
         eqy = 0.
         DO  j = 1, 10
           py = p(j)
           DO  i = 1, 10
             px = p(i)
             q1 = qs2val(px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,
     .                   dx,dy,rmax,rsq,a)
             CALL qs2grd(px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,
     .                   dx,dy,rmax,rsq,a,q,qx,qy,ier)
             IF (ier /= 0) GO TO 80
             IF (ABS(q1-q) > 3.*ABS(q)*eps) GO TO 90
             EQ = MAX(EQ,ABS(fq(px,py)-q))
             eqx = MAX(eqx,ABS(fx(px,py)-qx))
             eqy = MAX(eqy,ABS(fy(px,py)-qy))
             WRITE(0,*) 'TEST:',q,fq(px,py)
           END DO
         END DO

! PRINT ERRORS AND THE RATIO EQ/EPS.

         rq = EQ / eps
         WRITE (lout,5000)
         WRITE (lout,5100) EQ, rq
         WRITE (lout,5200) eqx
         WRITE (lout,5300) eqy
         STOP
       END IF

! ERROR IN QSHEP2

       WRITE (lout,5400) ier
       STOP

! ERROR IN QS2GRD

80        WRITE (lout,5500) ier
       STOP

! VALUES RETURNED BY QS2VAL AND QS2GRD DIFFER BY A RELATIVE
!   AMOUNT GREATER THAN 3*EPS.

90        WRITE (lout,5600) q1, q
       STOP

5000        FORMAT(//' ','MAXIMUM ABSOLUTE ERRORS IN THE INTERPOLANT '//
     .                   'Q AND PARTIAL'// 
     .      ' DERIVATIVES QX AND QY RELATIVE TO MACHINE PRECISION EPS'//  
     .      t12, 'FUNCTION   MAX ERROR   MAX ERROR/EPS'/)
5100        FORMAT (t15, 'Q       ', e10.3, '       ', f5.2)
5200        FORMAT (t15, 'QX      ', e10.3)
5300        FORMAT (t15, 'QY      ', e10.3)
5400        FORMAT (///' *** ERROR IN QSHEP2 -- IER =', i2, ' ***')
5500        FORMAT (///' *** ERROR IN QS2GRD -- IER =', i2, ' ***')
5600        FORMAT (///' *** ERROR -- INTERPOLATED VALUES ',  
     .       'Q1 (QS2VAL) AND Q2 (QS2GRD) DIFFER --',  
     .       t7, 'Q1 = ', e21.14, '     Q2 = ', e21.14)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE ClearTetArrays
      USE mod_geometry
      IMPLICIT none

      ngrp = 0
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
c
      SUBROUTINE GetTetVtx(newvtx,ivtx,status)
      USE mod_geometry
      IMPLICIT none

      REAL*8  newvtx(3)
      INTEGER ivtx,status

      status = 0
  
      IF (ivtx.LE.0.OR.ivtx.GT.nvtx) THEN
        status = -1
        RETURN
      ENDIF

      newvtx(1:3) = vtx(1:3,ivtx)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
      INTEGER FUNCTION GetTetVtxInd(iobj,iside,ivtx)
      USE mod_geometry
      IMPLICIT none

      INTEGER iobj,iside,ivtx

      INTEGER isrf

      isrf = ABS(obj(iobj)%iside(iside))

      GetTetVtxInd = srf(isrf)%ivtx(ivtx)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE HackTriangles(nhack)
      USE mod_eirene04
      IMPLICIT none
  
      INTEGER nhack

      ntri = nhack

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      REAL FUNCTION GetTetCentre(iobj)
      USE mod_geometry
      IMPLICIT none

      INTEGER iobj,isid,isrf,ivtx
      REAL*8  cen(3),count

      count = 0.0D0
      cen = 0.0D0
      DO isid = 1, obj(iobj)%nside
        isrf = ABS(obj(iobj)%iside(isid))
        DO ivtx = 1, srf(isrf)%nvtx
          count = count + 1.0D0
          cen(1:3) = cen(1:3) + vtx(1:3,srf(isrf)%ivtx(ivtx))
        ENDDO
      ENDDO
      cen = cen / count

      GetTetCentre = SNGL(cen(2))

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      INTEGER FUNCTION GetSurfaceIndex(iobj,iside)
      USE mod_geometry
      USE mod_sol28_global
      IMPLICIT none

      INTEGER iobj,iside,isrf

      isrf = ABS(obj(iobj)%iside(iside))
      isrf = srf(isrf)%index(IND_SURFACE)

      GetSurfaceIndex = isrf

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE GetNextTet(newobj,nsurface,ielement,option,status)
      USE mod_out985
      USE mod_geometry
      USE mod_eirene06_locals
      IMPLICIT none

      TYPE(type_3D_object) :: newobj
      INTEGER, INTENT(IN)  :: nsurface,ielement,option
      INTEGER, INTENT(OUT) :: status

      INTEGER GetSurfaceIndex
      REAL    GetTetCentre

      INTEGER ivolume,i1,i2,count,isrc,ivtx,iobj,iside,isrf,origin
      REAL, ALLOCATABLE :: ycen(:)

      SAVE

      status = 0

      IF (nsurface.EQ.-1) THEN
c...    Needed so that Dalpha could be loaded using LoadTriangleData
c       from AssignEmissionData... scary and will not doubt cause 
c       problems:
        WRITE(0,*)
        WRITE(0,*) '--------------------------------------------------'
        WRITE(0,*) '    HACKING ntri TO nobj SO LOADTRI... WORKS!     '
        WRITE(0,*) '--------------------------------------------------'
        WRITE(0,*)
        CALL HackTriangles(nobj)
        ivolume = 0
c        ALLOCATE(tdata(nobj))
c        CALL LoadTriangleData(6,1,7 ,1,tdata)  ! Dalpha
        ALLOCATE(ycen(nobj))
        DO i1 = 1, nobj
          ycen(i1) = GetTetCentre(i1)
        ENDDO
        iobj = 0
        RETURN
      ENDIF

c      iobj = iobj + 1
c      IF (iobj.GT.nobj) THEN
c        IF (ALLOCATED(tdata)) DEALLOCATE(tdata)
c        IF (ALLOCATED(ycen )) DEALLOCATE(ycen )
c        status = -1
c        RETURN
c      ENDIF

c      DO WHILE ((grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID).OR.
c     .          (ycen(iobj).LT.-0.90))
c     .          (ycen(iobj).LT.-0.50).OR.(ycen(iobj).GT.0.50))
c     .          (ycen(iobj).LT.-0.20).OR.(ycen(iobj).GT.0.20))
c     .          (tdata(iobj).LT.0.5E+23))
c     DO WHILE (obj(iobj)%segment(1).EQ.0)  ! *** HACK ***

c      DO WHILE (grp(obj(iobj)%group)%origin.NE.GRP_VACUUM_GRID)
cc      DO WHILE (grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID)
c        iobj = iobj + 1
c        IF (iobj.GT.nobj) THEN
c          IF (ALLOCATED(tdata)) DEALLOCATE(tdata)
c          IF (ALLOCATED(ycen )) DEALLOCATE(ycen )
c          status = -1
c          RETURN
c        ENDIF
c      ENDDO

c...  Select next valid tetrahedron:
      DO WHILE (.TRUE.) 
        iobj = iobj + 1
        IF (iobj.GT.nobj) THEN
c          IF (ALLOCATED(tdata)) DEALLOCATE(tdata)
          IF (ALLOCATED(ycen )) DEALLOCATE(ycen )
          status = -1
          RETURN
        ENDIF
        origin = grp(obj(iobj)%group)%origin
        IF (option.EQ.0.AND.origin.EQ.GRP_MAGNETIC_GRID) EXIT
        IF (option.EQ.1.AND.origin.EQ.GRP_VACUUM_GRID  ) EXIT
        IF (option.EQ.2                                ) EXIT
      ENDDO

      ivolume = ivolume + 1

      newobj%index        = ielement  ! nobj
      newobj%type         = OP_INTEGRATION_VOLUME
      newobj%subtype      = OP_EIRENE_GRID
      newobj%ivolume      = 1 ! ivolume  ! For line shape calculations, not sure what to do here...
      newobj%mode         = 0      
      newobj%surface      = 0        ! Transparent...?
      newobj%phi          = obj(iobj)%phi
      newobj%wedge1       = 0
      newobj%wedge2       = 0
      newobj%colour       = 3
      newobj%orientation  = 1      ! CW
      newobj%ik           = obj(iobj)%index(IND_IK)
      newobj%ir           = obj(iobj)%index(IND_IR)
      newobj%in           = iobj  
      newobj%nside        = 4
      newobj%iside(1,1:2) = nsurface + 1
      newobj%iside(2,1:2) = nsurface + 2
      newobj%iside(3,1:2) = nsurface + 3
      newobj%iside(4,1:2) = nsurface + 4
      newobj%gsur(1:4)    = GT_TD
      newobj%tsur(1:4)    = SP_GRID_SURFACE  ! Some to be replace with SP_GRID_BOUNDARY...
      newobj%reflec(1:4)  = 0
      newobj%nmap         = 1
      newobj%imap(1,1:4)  = obj(iobj)%omap(1:4)
      newobj%isur         = -1
      newobj%quantity(1)  = 1.0
c      IF (obj(iobj)%segment(1).EQ.0) newobj%quantity(1) = 0.0  ! *** HACK ***
c..   Defunct:  ! ...really?
      newobj%nsur         = 0
      newobj%ipts(2,1)    = 0

c...  Store the associated non-default standard surface index from 
c     the list in the EIRENE input file:
      DO iside = 1, obj(iobj)%nside
        newobj%esurf(iside) = GetSurfaceIndex(iobj,iside)
      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE Wrapper_LoadObjects(fname,status)
      USE mod_geometry
      IMPLICIT none

      INTEGER   status
      CHARACTER fname*(*)

      CALL LoadObjects(TRIM(fname),status)


c...  Dump tetrahedron data to an ASCII file:







      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE Wrapper_ClearObjects
      USE mod_geometry
      IMPLICIT none

      IF (ALLOCATED(obj)) DEALLOCATE(obj)
      IF (ALLOCATED(srf)) DEALLOCATE(srf)
      IF (ALLOCATED(vtx)) DEALLOCATE(vtx)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
      SUBROUTINE ProcessTetrahedronGrid(ielement)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

c...  Input:
      INTEGER ielement

      INTEGER AddVertex,AddSurface,GetTetVtxInd

      TYPE(type_3D_object) newobj
      TYPE(type_surface) newsrf
      REAL*8  newvtx(3)
      INTEGER status,i1,idum1

c...  For connection map:
      INTEGER nlist,i2,i3,i4,i5,i6,c,iobj,isid,isrf,isrf1,isrf2,
     .        iside,iside1,iobj1,ntet,nmatch,option
      INTEGER, ALLOCATABLE :: vsum(:,:),ilist(:), itet(:)
      REAL, ALLOCATABLE :: yobj(:)
      REAL    minphi,maxphi,phi,dphi,y,dy,miny,maxy,GetTetCentre


      WRITE(0,*) 'Loading tetrahedron grid',nobj
       

      CALL Wrapper_LoadObjects(TRIM(opt%obj_fname(ielement)),status)
c      CALL Wrapper_LoadObjects('tetrahedrons.raw',status)
      IF (status.NE.0) CALL ER('Wrapper_LoadObjects','Unable '//
     .                         'to find grid file',*99)



c.... Load all the vertices:
      status = 0
      i1 = 0
      DO WHILE(status.EQ.0)
        i1 = i1 + 1
        CALL GetTetVtx(newvtx,i1,status)
        IF (status.EQ.0) idum1 = AddVertex(newvtx)
      ENDDO

c...  Load all fluid grid tetrahedrons:
      option = opt%obj_option(ielement)
      status = 0
      CALL GetNextTet(newobj,-1,-1,-1,status)
      DO WHILE (status.EQ.0) 
        CALL GetNextTet(newobj,nsrf,ielement,option,status)

        IF (status.EQ.0) THEN
          DO i1 = 1, newobj%nside
            newsrf%type = SP_PLANAR_POLYGON
            newsrf%nvtx = 3       
            newsrf%ivtx(1) = GetTetVtxInd(newobj%in,i1,1)
            newsrf%ivtx(2) = GetTetVtxInd(newobj%in,i1,2)
            newsrf%ivtx(3) = GetTetVtxInd(newobj%in,i1,3)
            idum1 = AddSurface(newsrf)
          ENDDO
          IF (nobj+1.GT.MAX3D) 
     .      CALL ER('ProcessTetrahedronGrid','Insufficient array '//
     .              'bounds for all objects (A)',*99)    
          nobj = nobj + 1
          obj(nobj) = newobj
        ENDIF
      ENDDO 

      CALL Wrapper_ClearObjects

      WRITE(0,*) 'TETRAHEDRONS NSRF:',nsrf
      WRITE(0,*) 'TETRAHEDRONS NOBJ:',nobj
      WRITE(0,*) 'TETRAHEDRONS NVTX:',nvtx
c
c...  Need to build connection map, again:
      IF (.TRUE.) THEN
        WRITE(0,*) '  CALCULATING TETRAHEDRON CONNECTION MAP'
c...    Create map of original tetrahedron grid object index to the
c       new RAY object index:
        ntet = 0
        DO iobj = 1, nobj
          ntet = MAX(ntet,obj(iobj)%in)
        ENDDO
        ALLOCATE(itet(ntet))
        itet = -1
        DO iobj = 1, nobj
          itet(obj(iobj)%in) = iobj
        ENDDO
c...    Create (poor) surface discriminator:
        ALLOCATE(vsum(4,nobj))
        WRITE(0,*) '  CALCULATING VSUM'
        DO iobj = 1, nobj
          DO iside = 1, obj(iobj)%nside
            isrf = obj(iobj)%iside(iside,1)  ! Assuming only one surface per side...
            vsum(iside,iobj) = SUM(srf(isrf)%ivtx(1:3))
          ENDDO
        ENDDO
        WRITE(0,*) '  DONE'
c...    Build map:
        DO iobj = 1, nobj
          DO iside = 1, obj(iobj)%nside

            IF (obj(iobj)%imap(1,iside).NE.0   .AND.
     .          obj(iobj)%imap(1,iside).LE.ntet) THEN
              obj(iobj)%imap(1,iside) = itet(obj(iobj)%imap(1,iside))
            ELSE
              obj(iobj)%imap(1,iside) = -1
            ENDIF

            IF (obj(iobj)%imap(1,iside).EQ.-1) THEN
c...          No neighbouring tetrahedron:
              obj(iobj)%tsur(  iside) = SP_GRID_BOUNDARY
              obj(iobj)%imap(1,iside) = iobj
              obj(iobj)%isur(1,iside) = iside
            ELSE
c...          Find out which side the neigbouring tetrahedron that
c             the current side corresponds to:
              iobj1 = obj(iobj)%imap(1,iside)
              nmatch = 0
              DO iside1 = 1, obj(iobj1)%nside
                IF (vsum(iside,iobj).EQ.vsum(iside1,iobj1)) THEN
                  nmatch = nmatch + 1
                  obj(iobj)%isur(1,iside) = iside1                  
                ENDIF
              ENDDO
              IF (nmatch.EQ.0) 
     .          CALL ER('ProcessTetrahedronGrid','No side found '//
     .                  'when building connection map',*98)
              IF (nmatch.GT.1) 
     .          CALL ER('ProcessTetrahedronGrid','More than one '//
     .                  'side found, need to be more careful',*98)
            ENDIF
          ENDDO
        ENDDO
        DEALLOCATE(itet)
        DEALLOCATE(vsum)
        WRITE(0,*) '  DONE'
      ENDIF




c      DO iobj = 1, nobj
c        IF (obj(iobj)%ivolume.NE.2) CYCLE  
c        IF (obj(iobj)%ir     .NE.9) CYCLE  
c        WRITE(0,*) 'OBJ:',iobj,obj(iobj)%ik, 
c     .    obj(iobj)%tsur(1).EQ.SP_GRID_BOUNDARY
c      ENDDO
c      STOP 'sdfsdf'



      RETURN
 98   WRITE(0,*) ' FNAME= "'//TRIM(opt%obj_fname(ielement))//'"'
      WRITE(0,*) ' IOBJ,ISIDE = ',iobj,iside
      WRITE(0,*) ' NOBJ       = ',nobj
      WRITE(0,*) ' VSUM       = ',vsum(iside,iobj)
      WRITE(0,*) ' VSUM1      = ',vsum(1:4,obj(iobj)%imap(1,iside))
 99   WRITE(0,*) ' MAX3D      = ',MAX3D
      WRITE(0,*) ' NOBJ       = ',nobj
      STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE ProcessTriangleGrid(ielement)
      USE mod_out985
      USE mod_out985_variables
      USE mod_eirene06_parameters
      USE mod_eirene06
c      USE mod_eirene04
      use mod_params
      use mod_cgeom
      use mod_slcom
      IMPLICIT none

      INTEGER ielement

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'slcom'

      LOGICAL PointOnLine

      REAL*8     TOL        
      PARAMETER (TOL=1.0D-05)

      INTEGER v1,v2,itri,isector,nsector,ivol,i1,imap
      REAL    angle,dangle,dangle2,ang,rcen,zcen
      REAL*8  x1,z1,x2,z2,p1(3,8),p2(3,8),s1,s2,t1,t2,x(3),y(3)

c...  Load triangle .raw file:
      CALL LoadTriangles_06
c      CALL LoadTriangles

      nsector = opt%obj_nsector
      IF (nsector.EQ.-1) nsector = eirntorseg

c...  Just pulling wall surfaces:
c
c       opt%ob_trigrd = 3 - wall+target, toroidally discrete
c                       4 - wall
c                       5 - wall+target, toroidally continous
c                       6 - wall 
c                       7 - wall+target, toroidally discrete, non-symmetric
c                       8 - wall
c
      IF ((ielement.NE.0.AND.
     .     (opt%obj_option(ielement).EQ.3.OR.
     .      opt%obj_option(ielement).EQ.4.OR.
     .      opt%obj_option(ielement).EQ.5.OR.
     .      opt%obj_option(ielement).EQ.6.OR.
     .      opt%obj_option(ielement).EQ.7.OR.
     .      opt%obj_option(ielement).EQ.8))) THEN

        IF ((ielement.NE.0.AND.
     .       (opt%obj_option(ielement).EQ.5.OR.
     .        opt%obj_option(ielement).EQ.6))) THEN
c...      Toroidally continuous surfaces:
          nsector = 1
        ELSEIF ((ielement.NE.0.AND.
     .           (opt%obj_option(ielement).EQ.3.OR.
     .            opt%obj_option(ielement).EQ.4.OR.
     .            opt%obj_option(ielement).EQ.7.OR.
     .            opt%obj_option(ielement).EQ.8))) THEN
c...      Toroidally discrete:
          nsector = opt%obj_nsector
          IF (nsector.EQ.-1) nsector = eirntorseg
          dangle = 360.0 / REAL(nsector) / RADDEG
          nsector = NINT(REAL(nsector) * 
     .                  (opt%obj_angle_end-opt%obj_angle_start) / 360.0)
        ELSE
          CALL ER('ProcessTriangleGrid','Unknown triangle grid '//
     .            'option',*99)
        ENDIF

        dangle2 = dangle

        DO isector = 1, nsector

          dangle = dangle2  ! For ...

          ang = REAL(isector - 1) * dangle + opt%obj_angle_start/RADDEG

          IF ((ielement.NE.0.AND.
     .         (opt%obj_option(ielement).EQ.7.OR.
     .          opt%obj_option(ielement).EQ.8))) THEN
            IF (MOD(isector,2).EQ.1) THEN
              dangle = dangle2 * 1.60
            ELSE
              dangle = dangle2 * 0.40
            ENDIF
          ENDIF

          DO itri = 1, ntri
            DO v1 = 1, 3
              imap = tri(itri)%map(v1) 
              IF (.NOT.
     .            (tri(itri)%sur(v1).NE.0.AND.
c... 
     .             ((tri(itri)%map(v1).EQ.0            .AND.               ! Need surface type identifier...
     .               tri(itri)%type   .NE.MAGNETIC_GRID).OR.
c... 
     .              (
     .               ((ielement.NE.0.AND.
     .                 (opt%obj_option(ielement).EQ.3.OR.
     .                  opt%obj_option(ielement).EQ.5.OR.
     .                  opt%obj_option(ielement).EQ.7))).AND.
     .               tri(itri)%map(v1)        .EQ.0            .AND.
     .               tri(itri)%type           .EQ.MAGNETIC_GRID.AND.
     .               tri(itri)%sideindex(2,v1).NE.0).OR.                   
c... 
     .              (imap.NE.0.AND.
     .               tri(itri       )%type.NE.MAGNETIC_GRID.AND.
     .               tri(MAX(1,imap))%type.NE.MAGNETIC_GRID))))
     .          CYCLE

c...          Filter:
              IF (opt%obj_n(ielement,1).NE.-1.AND.
     .            opt%obj_n(ielement,2).NE.-1.AND.
     .         ((tri(itri)%sideindex(3,v1).NE.0.AND.
     .           (opt%obj_n(ielement,1).NE.0.OR.
     .            opt%obj_n(ielement,2).NE.0)).OR.
     .          (tri(itri)%sideindex(4,v1).NE.0.AND.
     .           (tri(itri)%sideindex(4,v1).LT.opt%obj_n(ielement,1).OR.
     .            tri(itri)%sideindex(4,v1).GT.opt%obj_n(ielement,2)))))
     .          CYCLE

              IF (nobj+1.GT.MAX3D) 
     .          CALL ER('ProcessTriangleGrid','Insufficient array '//
     .                  'bounds for all objects B',*98)     

              nobj = nobj + 1

              obj(nobj)%index       = ielement
              obj(nobj)%type        = OP_EMPTY
              obj(nobj)%mode        = 0      
              obj(nobj)%surface     = 1      ! SOLID
              obj(nobj)%phi         = ang
              obj(nobj)%wedge1      = 0
              obj(nobj)%wedge2      = 0
              obj(nobj)%colour      = 1
              obj(nobj)%orientation = 1      ! CW
              obj(nobj)%ik          = itri
              obj(nobj)%ir          = 0
              obj(nobj)%in          = 0
              obj(nobj)%ivolume     = 0
              obj(nobj)%nsur        = 1
              IF ((ielement.NE.0.AND.
     .             (opt%obj_option(ielement).EQ.5.OR.
     .              opt%obj_option(ielement).EQ.6))) THEN
                obj(nobj)%gsur(1:1) = GT_TC
                obj(nobj)%nver      = 2
                obj(nobj)%tsur(1)   = SP_VESSEL_WALL
                IF (ielement.NE.0) THEN
                  obj(nobj)%reflec(1) = opt%obj_reflec(ielement)
                ELSE
c                  obj(nobj)%reflec(1) = opt%ob_trigrd_reflec
                ENDIF
                obj(nobj)%npts(1)   = 2
                obj(nobj)%ipts(1,1) = 1
                obj(nobj)%ipts(2,1) = 2
              ELSE
                obj(nobj)%gsur(1:1) = GT_TD
                obj(nobj)%nver      = 4
                obj(nobj)%tsur(1)   = SP_VESSEL_WALL
                IF (ielement.NE.0) THEN
                  obj(nobj)%reflec(1) = opt%obj_reflec(ielement)
                ELSE
c                  obj(nobj)%reflec(1) = opt%ob_trigrd_reflec
                ENDIF
                obj(nobj)%npts(1)   = 4
                obj(nobj)%ipts(1,1) = 1
                obj(nobj)%ipts(2,1) = 2
                obj(nobj)%ipts(3,1) = 3
                obj(nobj)%ipts(4,1) = 4
              ENDIF
              obj(nobj)%nmap(1) = 0
 
c...          Vertices:
              v2 = v1 + 1
              IF (v2.EQ.4) v2 = 1 

              IF ((ielement.NE.0.AND.
     .             (opt%obj_option(ielement).EQ.5.OR.
     .              opt%obj_option(ielement).EQ.6))) THEN
                obj(nobj)%v(1,1) = DBLE(ver(tri(itri)%ver(v1),1))
                obj(nobj)%v(2,1) = DBLE(ver(tri(itri)%ver(v1),2))
                obj(nobj)%v(3,1) = 0.0D0
                obj(nobj)%v(1,2) = DBLE(ver(tri(itri)%ver(v2),1))
                obj(nobj)%v(2,2) = DBLE(ver(tri(itri)%ver(v2),2))
                obj(nobj)%v(3,2) = 0.0D0
              ELSE
                obj(nobj)%v(1,1) = DBLE(ver(tri(itri)%ver(v1),1)) *   
     .                             DCOS(DBLE(-0.5*dangle))           
                obj(nobj)%v(2,1) = DBLE(ver(tri(itri)%ver(v1),2))   
                obj(nobj)%v(3,1) = DBLE(ver(tri(itri)%ver(v1),1)) * 
     .                             DSIN(DBLE(-0.5*dangle))           
                obj(nobj)%v(1,2) = DBLE(ver(tri(itri)%ver(v2),1)) *   
     .                             DCOS(DBLE(-0.5*dangle))           
                obj(nobj)%v(2,2) = DBLE(ver(tri(itri)%ver(v2),2))   
                obj(nobj)%v(3,2) = DBLE(ver(tri(itri)%ver(v2),1)) * 
     .                             DSIN(DBLE(-0.5*dangle))           
                obj(nobj)%v(1,3) = DBLE(ver(tri(itri)%ver(v2),1)) *   
     .                             DCOS(DBLE(+0.5*dangle))           
                obj(nobj)%v(2,3) = DBLE(ver(tri(itri)%ver(v2),2))   
                obj(nobj)%v(3,3) = DBLE(ver(tri(itri)%ver(v2),1)) * 
     .                             DSIN(DBLE(+0.5*dangle))           
                obj(nobj)%v(1,4) = DBLE(ver(tri(itri)%ver(v1),1)) *   
     .                             DCOS(DBLE(+0.5*dangle))           
                obj(nobj)%v(2,4) = DBLE(ver(tri(itri)%ver(v1),2))   
                obj(nobj)%v(3,4) = DBLE(ver(tri(itri)%ver(v1),1)) * 
     .                             DSIN(DBLE(+0.5*dangle))           
c...            Rotate vertices:
                DO i1 = 1, 4
                  x1 = obj(nobj)%v(1,i1)
                  z1 = obj(nobj)%v(3,i1)
                  obj(nobj)%v(1,i1) = DCOS(DBLE(ang)) * x1 -
     .                                DSIN(DBLE(ang)) * z1
                  obj(nobj)%v(3,i1) = DSIN(DBLE(ang)) * x1 +
     .                                DCOS(DBLE(ang)) * z1
                ENDDO
              ENDIF

            ENDDO
          ENDDO
        ENDDO

c...    Horizontal test surface for checking distortion: 
        IF (.FALSE.) THEN
          nobj = nobj + 1
          obj(nobj)%index       = ielement
          obj(nobj)%type        = OP_EMPTY
          obj(nobj)%mode        = 0      
          obj(nobj)%surface     = 1      ! SOLID
          obj(nobj)%wedge1      = 0
          obj(nobj)%wedge2      = 0
          obj(nobj)%colour      = 1
          obj(nobj)%orientation = 1      ! CW
          obj(nobj)%ik          = 0
          obj(nobj)%ir          = 0
          obj(nobj)%in          = ntri+1
          obj(nobj)%ivolume     = 0
          obj(nobj)%nsur        = 1
          obj(nobj)%gsur(1:1)   = GT_TD
          obj(nobj)%nver        = 2
          obj(nobj)%tsur(1)     = SP_VESSEL_WALL
          obj(nobj)%reflec(1)   = 0 
          obj(nobj)%npts(1)     = 4
          obj(nobj)%ipts(1,1)   = 1
          obj(nobj)%ipts(2,1)   = 2
          obj(nobj)%ipts(3,1)   = 3
          obj(nobj)%ipts(4,1)   = 4
          obj(nobj)%nmap(1) = 0
c...      Vertices:
          obj(nobj)%v(1,1) =  0.5D0
          obj(nobj)%v(2,1) =  1.3D0
          obj(nobj)%v(3,1) =  1.0D0
          obj(nobj)%v(1,2) =  0.5D0
          obj(nobj)%v(2,2) =  1.3D0
          obj(nobj)%v(3,2) = -1.0D0
          obj(nobj)%v(1,3) =  0.5D0
          obj(nobj)%v(2,3) =  1.2D0
          obj(nobj)%v(3,3) = -1.0D0
          obj(nobj)%v(1,4) =  0.5D0
          obj(nobj)%v(2,4) =  1.2D0
          obj(nobj)%v(3,4) =  1.0D0
        ENDIF
c     ------------------------------------------------------------------
      ELSEIF (opt%obj_option(ielement).EQ.2.OR.
     .        opt%obj_option(ielement).EQ.9) THEN
        ivol = 0
        DO itri = 1, ntri
          IF (tri(itri)%type.NE.MAGNETIC_GRID) CYCLE

          IF (nobj+1.GT.MAX3D) 
     .      CALL ER('BuildObjects','Insufficient array bounds '//
     .              'for all objects C',*98)     
          nobj = nobj + 1
          obj(nobj)%index       = ielement ! nobj

          ivol = ivol + 1
          obj(nobj)%type    = OP_INTEGRATION_VOLUME
          obj(nobj)%subtype = OP_EIRENE_GRID
          obj(nobj)%ivolume = ivol

          obj(nobj)%mode        = 0      
          obj(nobj)%surface     = 1      ! SOLID
          obj(nobj)%wedge1      = 0
          obj(nobj)%wedge2      = 0
          obj(nobj)%colour      = 3
          obj(nobj)%orientation = 1      ! CW
          obj(nobj)%ik          = tri(itri)%index(1)
          obj(nobj)%ir          = tri(itri)%index(2)
          obj(nobj)%in          = itri
          obj(nobj)%nsur        = 3  
          obj(nobj)%gsur(1:5)   = GT_TC
          obj(nobj)%nver        = 3
          obj(nobj)%quantity    = 1.0
c...    
          DO v1 = 1, 3
            p1(1,v1  ) = DBLE(ver(tri(itri)%ver(v1),1))   ! *** NEED TO CHANGE VER(*,1:3) to VER(1:3,*) ***
            p1(2,v1  ) = DBLE(ver(tri(itri)%ver(v1),2))   
            p1(3,v1  ) = 0.0D0
c            p1(1,v1+3) = DBLE(ver(tri(itri)%ver(v1),1))  
c            p1(2,v1+3) = DBLE(ver(tri(itri)%ver(v1),2))  
c            p1(3,v1+3) = 0.0D0
          ENDDO
c...      Assign vertices to object:
c          IF (obj(nobj)%ik.EQ.66.AND.obj(nobj)%ir.EQ.56) THEN
c            WRITE(0,*) 'P1:',p1(1,1:3)
c            WRITE(0,*) '  :',p1(2,1:3)
c            STOP 'fsdf'
c          ENDIF

          DO v1 = 1, obj(nobj)%nver
            obj(nobj)%v(1:3,v1) = p1(1:3,v1)
          ENDDO
c...
          DO v1 = 1, 3

            IF (tri(itri)%sur(v1).NE.0) THEN        
c...          Triangle surface is on a surface (magnetic or vessel wall):
              imap = tri(itri)%map(v1)
              IF (opt%obj_option(ielement).EQ.2.AND.
     .            tri(itri)%map  (v1).EQ.0    .AND.
     .            tri(itri)%index(2 ).GE.irsep) THEN
c...            Vessel wall surface:
                obj(nobj)%tsur(v1) = SP_VESSEL_WALL  
                obj(nobj)%reflec(v1) = opt%obj_reflec(ielement)
                obj(nobj)%nmap(v1) = 1
                obj(nobj)%imap(1,v1) = nobj
                obj(nobj)%isur(1,v1) = v1 ! 2  ! *** should this really be a 2 for some reason? ***
                obj(nobj)%rsur(3,v1) = tri(itri)%sideindex(3,v1)   ! xVESM wall index
                obj(nobj)%rsur(4,v1) = tri(itri)%sideindex(4,v1)   ! Additional surface index
              ELSE
c...            Magnetic surface (grid boundary):
                obj(nobj)%tsur(v1) = SP_GRID_BOUNDARY 
                obj(nobj)%reflec(v1) = 0
                obj(nobj)%nmap(v1) = 1
                obj(nobj)%imap(1,v1) = nobj
                obj(nobj)%isur(1,v1) = v1    ! Surface maps to itself... wise...?
              ENDIF
            ELSE
c...          Triangle mesh boundary surface only:
              obj(nobj)%tsur(v1) = SP_GRID_SURFACE
              obj(nobj)%reflec(v1) = 0
              obj(nobj)%nmap(v1) = 1
              obj(nobj)%imap(1,v1) = tri(itri)%map(v1)
              obj(nobj)%isur(1,v1) = tri(itri)%sid(v1)
            ENDIF

            IF     (v1.EQ.1) THEN
              obj(nobj)%npts(v1) = 2
              obj(nobj)%ipts(1,v1) = 1
              obj(nobj)%ipts(2,v1) = 2
            ELSEIF (v1.EQ.2) THEN
              obj(nobj)%npts(v1) = 2
              obj(nobj)%ipts(1,v1) = 2
              obj(nobj)%ipts(2,v1) = 3
            ELSEIF (v1.EQ.3) THEN
              obj(nobj)%npts(v1) = 2
              obj(nobj)%ipts(1,v1) = 3
              obj(nobj)%ipts(2,v1) = 1
            ENDIF
          ENDDO


c          IF (obj(nobj)%ik.EQ.66.AND.obj(nobj)%ir.EQ.56) THEN
c            WRITE(0,*) 'P:',p1(1,1:3)
c            WRITE(0,*) ' :',p1(2,1:3)
c            DO v1 = 1, 3
c              WRITE(0,*) obj(nobj)%v(1,obj(nobj)%ipts(1,v1))
c              WRITE(0,*) obj(nobj)%v(1,obj(nobj)%ipts(2,v1))
c              WRITE(0,*) obj(nobj)%v(2,obj(nobj)%ipts(1,v1))
c              WRITE(0,*) obj(nobj)%v(2,obj(nobj)%ipts(2,v1))
c            ENDDO
c          ENDIF

        ENDDO
c     ------------------------------------------------------------------
      ELSE
        dangle = 360.0 / REAL(nsector) / RADDEG
        DO isector = 1, nsector
          ang = REAL(isector - 1) * dangle
          ivol = 0
          DO itri = 1, ntri
            IF (nobj+1.GT.MAX3D) 
     .        CALL ER('BuildObjects','Insufficient array bounds '//
     .                'for all objects D',*98)     

            nobj = nobj + 1

            obj(nobj)%index       = ielement  ! nobj
            IF (tri(itri)%type.EQ.MAGNETIC_GRID) THEN
              ivol = ivol + 1
              obj(nobj)%type    = OP_INTEGRATION_VOLUME
              obj(nobj)%subtype = OP_EIRENE_GRID
              obj(nobj)%ivolume = ivol
            ELSE
              obj(nobj)%type = OP_EMPTY
              obj(nobj)%ivolume = 0
            ENDIF
            obj(nobj)%mode        = 0      
            obj(nobj)%surface     = 1      ! SOLID
            obj(nobj)%phi         = ang
            obj(nobj)%wedge1      = 0
            obj(nobj)%wedge2      = 0
            obj(nobj)%colour      = 3
            obj(nobj)%orientation = 1      ! CW
            obj(nobj)%ik          = itri
            obj(nobj)%ir          = 0
            obj(nobj)%in          = 0
            obj(nobj)%nsur        = 5  ! 6 bug
            IF (opt%obj_option(ielement).EQ.2) THEN
              obj(nobj)%gsur(1:5) = GT_TC
              obj(nobj)%nver      = 8
            ELSE
              obj(nobj)%gsur(1:5) = GT_TD
              obj(nobj)%nver      = 8
            ENDIF


            obj(nobj)%quantity = 1.0
c            IF (obj(nobj)%type.EQ.OP_INTEGRATION_VOLUME) THEN
c              obj(nobj)%quantity(1) = 1.0    ! *TEMP*
c            ELSE
c              obj(nobj)%quantity(1) = 0.0
c            ENDIF

c            IF (tri(itri)%type.EQ.MAGNETIC_GRID.AND.
c     .          tri(itri)%index(2).EQ.8) THEN
c              obj(nobj)%quantity(2) = 1.0    
c            ELSE
c              obj(nobj)%quantity(2) = 0.0    
c            ENDIF

c...      
            DO v1 = 1, 3
              p1(1,v1  ) = DBLE(ver(tri(itri)%ver(v1),1))   ! *** NEED TO CHANGE VER(*,1:3) to VER(1:3,*) ***
              p1(2,v1  ) = DBLE(ver(tri(itri)%ver(v1),2))   ! *** NEED TO CHANGE VER(*,1:3) to VER(1:3,*) ***
              p1(3,v1  ) = DBLE(ver(tri(itri)%ver(v1),1)) * ! *** NOT USING Z-COORDINATE DATA *** 
     .                     DTAN(DBLE(-0.5*dangle))           
     .                 
              p1(1,v1+3) = DBLE(ver(tri(itri)%ver(v1),1))  
              p1(2,v1+3) = DBLE(ver(tri(itri)%ver(v1),2))  
              p1(3,v1+3) = DBLE(ver(tri(itri)%ver(v1),1)) * 
     .                     DTAN(DBLE(+0.5*dangle))           
            ENDDO
c...        Assign vertices to object:
            DO v1 = 1, 6
              obj(nobj)%v(1:3,v1) = p1(1:3,v1)
            ENDDO
c...        Rotate vertices:   ! *** WRONG WAY? ***
            DO v1 = 1, 6
              x1 = obj(nobj)%v(1,v1)
              z1 = obj(nobj)%v(3,v1)
              obj(nobj)%v(1,v1) = DCOS(DBLE(ang)) * x1 -
     .                            DSIN(DBLE(ang)) * z1
              obj(nobj)%v(3,v1) = DSIN(DBLE(ang)) * x1 +
     .                            DCOS(DBLE(ang)) * z1
            ENDDO
c...

            IF (isector.EQ.1) THEN   ! *** CHECK IF NOT FULL TORUS ***
              obj(nobj)%tsur(1) = SP_GRID_SURFACE
              obj(nobj)%reflec(1) = 0
              obj(nobj)%nmap(1) = 1
              obj(nobj)%imap(1,1) = nobj + ntri * (nsector - 1)
              obj(nobj)%isur(1,1) = 5
            ELSE
              obj(nobj)%tsur(1) = SP_GRID_SURFACE
              obj(nobj)%reflec(1) = 0
              obj(nobj)%nmap(1) = 1
              obj(nobj)%imap(1,1) = nobj - ntri
              obj(nobj)%isur(1,1) = 5
            ENDIF
            obj(nobj)%npts(1) = 3
            obj(nobj)%ipts(1,1) = 1
            obj(nobj)%ipts(2,1) = 2
            obj(nobj)%ipts(3,1) = 3

            DO v1 = 1, 3

c... Oy, need to sort this out.  Make sure %map is 0 for triangles on opaque surfaces... and
c    include some kind of identifier for the surface that the triangle side maps onto ---
c    perhaps the index of the surface in the xVESM and ADD SUR lists, along with the type of 
c    surface that it is... 

              IF (tri(itri)%sur(v1).NE.0) THEN        
c...            Triangle surface is on a surface (magnetic or vessel wall):
                IF ((tri(itri)%map(v1).EQ.0            .AND.               ! Need surface type identifier...
     .               tri(itri)%type   .NE.MAGNETIC_GRID).OR.
     .              (tri(itri)%map(v1)        .EQ.0            .AND.
     .               tri(itri)%type           .EQ.MAGNETIC_GRID.AND.
     .               tri(itri)%sideindex(2,v1).NE.0).OR.                   ! Target check
     .              (tri(itri)%map(v1)          .NE.0            .AND.
     .               tri(itri             )%type.NE.MAGNETIC_GRID.AND.
     .               tri(tri(itri)%map(v1))%type.NE.MAGNETIC_GRID)) THEN
c...              Vessel wall surface:
                  obj(nobj)%tsur(v1+1) = SP_VESSEL_WALL  
                 IF (ielement.NE.0) THEN
                    obj(nobj)%reflec(v1+1) = opt%obj_reflec(ielement)
                  ELSE
C                    obj(nobj)%reflec(v1+1) = opt%ob_trigrd_reflec
                  ENDIF
c                  obj(nobj)%reflec(v1+1) = opt%ob_trigrd_reflec
                  obj(nobj)%nmap(v1+1) = 1
                  obj(nobj)%imap(1,v1+1) = nobj
                  obj(nobj)%isur(1,v1+1) = 2  ! *** should this really be a 2 for some reason? ***
                  obj(nobj)%rsur(3,v1+1) = tri(itri)%sideindex(3,v1)   ! xVESM wall index
                  obj(nobj)%rsur(4,v1+1) = tri(itri)%sideindex(4,v1)   ! Additional surface index
                ELSE
c...              Magnetic surface (grid boundary):
                  obj(nobj)%tsur(v1+1) = SP_GRID_BOUNDARY  ! *** TRUE? ***
                  IF (ielement.NE.0) THEN
                    obj(nobj)%reflec(v1+1) = opt%obj_reflec(ielement)
                  ELSE
C                    obj(nobj)%reflec(v1+1) = opt%ob_trigrd_reflec
                  ENDIF
C                  obj(nobj)%reflec(v1+1) = opt%ob_trigrd_reflec
                  obj(nobj)%nmap(v1+1) = 1
                  obj(nobj)%imap(1,v1+1) = nobj
                  obj(nobj)%isur(1,v1+1) = v1 + 1  ! Surface maps to itself... wise...?
                ENDIF
              ELSE
c...           Triangle mesh boundary surface only:
               obj(nobj)%tsur(v1+1) = SP_GRID_SURFACE
               obj(nobj)%reflec(v1+1) = 0
               obj(nobj)%nmap(v1+1) = 1
               obj(nobj)%imap(1,v1+1)=tri(itri)%map(v1)+(isector-1)*ntri ! Untested
               obj(nobj)%isur(1,v1+1)=tri(itri)%sid(v1) + 1
              ENDIF

              IF     (v1.EQ.1) THEN
                obj(nobj)%npts(v1+1) = 4
                obj(nobj)%ipts(1,v1+1) = 2
                obj(nobj)%ipts(2,v1+1) = 5
                obj(nobj)%ipts(3,v1+1) = 4
                obj(nobj)%ipts(4,v1+1) = 1
c                obj(nobj)%ipts(1,v1) = 2
c                obj(nobj)%ipts(2,v1) = 1
c                obj(nobj)%ipts(3,v1) = 4
c                obj(nobj)%ipts(4,v1) = 5
              ELSEIF (v1.EQ.2) THEN
                obj(nobj)%npts(v1+1) = 4
                obj(nobj)%ipts(1,v1+1) = 3
                obj(nobj)%ipts(2,v1+1) = 6
                obj(nobj)%ipts(3,v1+1) = 5
                obj(nobj)%ipts(4,v1+1) = 2
c                obj(nobj)%ipts(1,v1) = 3
c                obj(nobj)%ipts(2,v1) = 2
c                obj(nobj)%ipts(3,v1) = 5
c                obj(nobj)%ipts(4,v1) = 6
              ELSEIF (v1.EQ.3) THEN
                obj(nobj)%npts(v1+1) = 4
                obj(nobj)%ipts(1,v1+1) = 1
                obj(nobj)%ipts(2,v1+1) = 4
                obj(nobj)%ipts(3,v1+1) = 6
                obj(nobj)%ipts(4,v1+1) = 3
c                 obj(nobj)%ipts(1,v1) = 1
c                 obj(nobj)%ipts(2,v1) = 3
c                 obj(nobj)%ipts(3,v1) = 6
c                 obj(nobj)%ipts(4,v1) = 4
              ENDIF

            ENDDO

            IF (isector.EQ.nsector) THEN ! *** CHECK IF NOT FULL TORUS! BIG SAVINGS?
              obj(nobj)%tsur(5) = SP_GRID_SURFACE
              obj(nobj)%reflec(5) = 0
              obj(nobj)%nmap(5) = 1
              obj(nobj)%imap(1,5) = nobj - ntri * (nsector - 1)
              obj(nobj)%isur(1,5) = 1
            ELSE
              obj(nobj)%tsur(5) = SP_GRID_SURFACE
              obj(nobj)%reflec(v1) = 0
              obj(nobj)%nmap(5) = 1
              obj(nobj)%imap(1,5) = nobj + ntri
              obj(nobj)%isur(1,5) = 1
            ENDIF
            obj(nobj)%npts(5) = 3
            obj(nobj)%ipts(1,5) = 6
            obj(nobj)%ipts(2,5) = 5
            obj(nobj)%ipts(3,5) = 4

          ENDDO

        ENDDO

      ENDIF

c...  Clear arrays:
      CALL DEALLOC_VERTEX  
      CALL DEALLOC_SURFACE   
      CALL DEALLOC_TRIANGLE

 98   RETURN
 99   STOP
      END
c
c ====================================================================== 
c
      SUBROUTINE ProcessMagneticGrid(ielement)
      USE mod_out985
      USE mod_out985_variables
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_dynam2
      use mod_slcom
      IMPLICIT none

      INTEGER ielement

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'pindata'
c     include 'dynam2'
c     INCLUDE 'slcom'


      LOGICAL PointOnLine

      REAL*8     TOL
      PARAMETER (TOL=1.0D-05)

      INTEGER i1,i2,i3,i4,ir,ik,id,ncells,isector,nsector,ike,ik1,ir1,
     .        ike1,id1,ivol,cntcell,maxmap
      LOGICAL check1,check2
      REAL    angle,dangle,ang,dangle2
      REAL*8  x1,z1,x2,z2,p1(3,8),p2(3,8),s1,s2,t1,t2,x(3),y(3)

      INTEGER, ALLOCATABLE :: ikinmap (:,:,:),irinmap (:,:,:),
     .                        ikoutmap(:,:,:),iroutmap(:,:,:),
     .                        ninmap(:,:),noutmap(:,:),indmap(:,:)

      INTEGER, PARAMETER :: MAXNMAP = 50

c...    Number of objects per toroidal segment:

c      WRITE(0,*) 'BUILDING RADIAL MAP'

      ALLOCATE(ninmap(MAXNKS,MAXNRS))
      ALLOCATE(ikinmap(MAXNKS,MAXNRS,MAXNMAP))
      ALLOCATE(irinmap(MAXNKS,MAXNRS,MAXNMAP))
      ALLOCATE(noutmap(MAXNKS,MAXNRS))
      ALLOCATE(ikoutmap(MAXNKS,MAXNRS,MAXNMAP))
      ALLOCATE(iroutmap(MAXNKS,MAXNRS,MAXNMAP))
      ALLOCATE(indmap(MAXNKS,MAXNRS))

      ninmap  = 0
      noutmap = 0
      

      ncells = 0

      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1
        DO ik = 1, ike

c...    Filter:
c          IF (zs(ik,ir).LT.1.0) CYCLE  
c          IF (zs(ik,ir).LT.-0.5.OR.zs(ik,ir).GT.0.5) CYCLE  
c          IF (zs(ik,ir).GT.-1.0) CYCLE
c          IF (zs(ik,ir).GT.-0.25) CYCLE

          ncells = ncells + 1

          id = korpg(ik,ir)
          x(1) = DBLE(rvertp(1,id))
          y(1) = DBLE(zvertp(1,id))
          x(2) = DBLE(rvertp(4,id))
          y(2) = DBLE(zvertp(4,id))

          ir1 = irins(ik,ir)
          ninmap(ik,ir) = 0
          IF (idring(ir1).NE.BOUNDARY) THEN
            ike1 = nks(ir1)
            IF (ir1.LT.irsep) ike1 = ike1 - 1
            DO ik1 = 1, ike1
              id1 = korpg(ik1,ir1)
              x(3) = DBLE(rvertp(2,id1))
              y(3) = DBLE(zvertp(2,id1))
              check1 = PointOnLine(x,y,s1,t1,4,.FALSE.) 
              x(3) = DBLE(rvertp(3,id1))
              y(3) = DBLE(zvertp(3,id1))
              check2 = PointOnLine(x,y,s2,t2,4,.FALSE.) 
              IF ((check1.AND.s1.GT.0.0D0+TOL.AND.s1.LT.1.0D0-TOL).OR.
     .            (check2.AND.s2.GT.0.0D0+TOL.AND.s2.LT.1.0D0-TOL).OR.
     .            (.NOT.check1.AND.check2.AND.s2.GT.1.0D0-TOL).OR.       ! Redundant? 
     .            (.NOT.check2.AND.check1.AND.s1.LT.0.0D0+TOL).OR.
     .            (check1.AND.check2)) THEN
                ninmap(ik,ir) = ninmap(ik,ir) + 1
                IF (ninmap(ik,ir).GT.MAXNMAP) 
     .            CALL ER('BuildObjects','Bounds error, '//  ! Note, this error can also come from not
     .                    'increase MAXNMAP A',*99)           ! loading the supplimental .raw files first
                ikinmap(ik,ir,ninmap(ik,ir)) = ik1
                irinmap(ik,ir,ninmap(ik,ir)) = ir1
              ENDIF
            ENDDO
          ENDIF

          x(1) = DBLE(rvertp(2,id))
          y(1) = DBLE(zvertp(2,id))
          x(2) = DBLE(rvertp(3,id))
          y(2) = DBLE(zvertp(3,id))
          ir1 = irouts(ik,ir)
          noutmap(ik,ir) = 0
          IF (idring(ir1).NE.BOUNDARY) THEN
            ike1 = nks(ir1)
            IF (ir1.LT.irsep) ike1 = ike1 - 1
            DO ik1 = 1, ike1
              id1 = korpg(ik1,ir1)
              x(3) = DBLE(rvertp(1,id1))
              y(3) = DBLE(zvertp(1,id1))
              check1 = PointOnLine(x,y,s1,t1,4,.FALSE.) 
              x(3) = DBLE(rvertp(4,id1))
              y(3) = DBLE(zvertp(4,id1))
              check2 = PointOnLine(x,y,s2,t2,4,.FALSE.) 
              IF ((check1.AND.s1.GT.0.0D0+TOL.AND.s1.LT.1.0D0-TOL).OR.
     .            (check2.AND.s2.GT.0.0D0+TOL.AND.s2.LT.1.0D0-TOL).OR.
     .            (.NOT.check1.AND.check2.AND.s2.GT.1.0D0-TOL).OR.
     .            (.NOT.check2.AND.check1.AND.s1.LT.0.0D0+TOL).OR.
     .            (check1.AND.check2)) THEN
                noutmap(ik,ir) = noutmap(ik,ir) + 1
                IF (noutmap(ik,ir).GT.MAXNMAP) 
     .            CALL ER('BuildObjects','Bounds error, '//
     .                    'increase MAXNMAP B',*99)
                ikoutmap(ik,ir,noutmap(ik,ir)) = ik1
                iroutmap(ik,ir,noutmap(ik,ir)) = ir1
              ENDIF
            ENDDO
          ENDIF

        ENDDO
      ENDDO

c *TEMP*
      maxmap = 0
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        DO ik = 1, nks(ir)

          maxmap = MAX(maxmap,MAX(ninmap(ik,ir),noutmap(ik,ir)))
c          IF(ninmap(ik,ir)

          IF (ir.LT.irsep.AND.ik.EQ.nks(ir)) CYCLE

c Filter
c          IF (zs(ik,ir).LT.1.0) CYCLE  
c          IF (zs(ik,ir).LT.-0.5.OR.zs(ik,ir).GT.0.5) CYCLE  
c          IF (zs(ik,ir).GT.-1.0) CYCLE
c          IF (zs(ik,ir).GT.-0.25) CYCLE

          DO i1 = 1, ninmap(ik,ir)
            IF (ikinmap(ik,ir,i1).EQ.0.OR.
     .          irinmap(ik,ir,i1).EQ.0) STOP 'sdfsdf A'
          ENDDO
          DO i1 = 1, noutmap(ik,ir)
            IF (ikoutmap(ik,ir,i1).EQ.0.OR.
     .          iroutmap(ik,ir,i1).EQ.0) STOP 'sdfsdf B'
          ENDDO
 
        ENDDO
      ENDDO
      
      WRITE(0,*) 'MAXMAP:',maxmap

c      STOP 'dasdsad'
c      WRITE(0,*) 'DONE'

      IF     ((ielement.NE.0.AND.
     .         (opt%obj_option(ielement).EQ.2.OR.
     .          opt%obj_option(ielement).EQ.4))) THEN
C     .        (opt%ob_stdgrd.EQ.2)) THEN
c...    Toroidally continuous surfaces:
        nsector = 1
        dangle = 0.0
      ELSEIF ((ielement.NE.0.AND.
     .         (opt%obj_option(ielement).EQ.1.OR.
     .          opt%obj_option(ielement).EQ.5.OR.
     .          opt%obj_option(ielement).EQ.3))) THEN
C     .        (opt%ob_stdgrd.EQ.1.OR.opt%ob_stdgrd.EQ.3)) THEN
        nsector = opt%obj_nsector

c        STOP 'CUSTOM TOROIDAL SECTION SPECIFIED'
        nsector = 12
        opt%obj_angle_start = 0.0
        opt%obj_angle_end = 3.0 * 360.0 / 48.0
!        opt%obj_angle_end = 360.0

        IF (nsector.EQ.-1) nsector = eirntorseg
        dangle = 360.0 / REAL(nsector) / RADDEG
        nsector = NINT(REAL(nsector) * 
     .                 (opt%obj_angle_end-opt%obj_angle_start) / 360.0)
      ELSE
        WRITE(0,*) 'OPTION: ',ielement,opt%obj_option(ielement)
        CALL ER('ProcessMagneticGrid','Unknown standard grid '//
     .          'option',*99)
      ENDIF

      dangle2 = dangle

      DO isector = 1, nsector

        dangle = dangle2  ! For opt%ob_stdgrd.EQ.3, temp...

        ang = REAL(isector - 1) * dangle + opt%obj_angle_start / RADDEG

        IF ((ielement.NE.0.AND.
     .       opt%obj_option(ielement).EQ.3)) THEN
c     .      opt%ob_stdgrd.EQ.3) THEN
          IF (MOD(isector,2).EQ.1) THEN
            dangle = dangle2 * 1.60
          ELSE
            dangle = dangle2 * 0.40
          ENDIF
        ENDIF

        ivol = 0

        DO ir = 1, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          ike = nks(ir)
          IF (ir.LT.irsep) ike = ike - 1
c...      Count the number of cells on this magnetic grid ring that will
c         make it onto the integration mesh, which isn't always equal to
c         IKE due to mesh filter/mask rules: 
          cntcell = 0
          DO ik = 1, ike
c Filter
c            IF (zs(ik,ir).LT.1.0) CYCLE  
c            IF (zs(ik,ir).LT.-0.5.OR.zs(ik,ir).GT.0.5) CYCLE  
c            IF (zs(ik,ir).GT.-1.0) CYCLE
c            IF (zs(ik,ir).GT.-0.25) CYCLE
            cntcell = cntcell + 1        
          ENDDO                          
c...
          DO ik = 1, ike

            IF (nobj+1.GT.MAX3D) 
     .        CALL ER('BuildObjects','Insufficient array bounds '//
     .                'for all objects E',*98)     

c Filter
c            IF (zs(ik,ir).LT.1.0) CYCLE  
c            IF (zs(ik,ir).LT.-0.5.OR.zs(ik,ir).GT.0.5) CYCLE   
c            IF (zs(ik,ir).GT.-1.0) CYCLE  
c            IF (zs(ik,ir).GT.-0.25) CYCLE

            ivol = ivol + 1

            nobj = nobj + 1

            obj(nobj)%index       = ielement  ! nobj
            obj(nobj)%type        = OP_INTEGRATION_VOLUME
            obj(nobj)%subtype     = OP_FLUID_GRID
            obj(nobj)%mode        = 0      
            obj(nobj)%surface     = 1      ! SOLID
            obj(nobj)%phi         = ang
            obj(nobj)%wedge1      = 0
            obj(nobj)%wedge2      = 0
            obj(nobj)%colour      = 3
            obj(nobj)%orientation = 1      ! CW
            obj(nobj)%ik          = ik
            obj(nobj)%ir          = ir
            obj(nobj)%in          = 0
            obj(nobj)%ivolume     = ivol

           IF ((ielement.NE.0.AND.
     .          (opt%obj_option(ielement).EQ.2.OR.
     .           opt%obj_option(ielement).EQ.4))) THEN
c     .         opt%ob_stdgrd.EQ.2) THEN
              obj(nobj)%nsur        = 4
              obj(nobj)%gsur(1:4)   = GT_TC
              obj(nobj)%nver        = 4
            ELSE
              obj(nobj)%nsur        = 6
              obj(nobj)%gsur(1:6)   = GT_TD
              obj(nobj)%nver        = 8
            ENDIF

            IF (.TRUE.) THEN
 
              obj(nobj)%quantity = 1.0

            ELSE
c              obj(nobj)%quantity(1) = 1.0
              obj(nobj)%quantity(1) = pinline(ik,ir,6,H_BGAMMA)
c              obj(nobj)%quantity(1) = pinalpha(ik,ir)
              obj(nobj)%quantity(2) = 10.0D0+
     .                           LOG10(MAX(1.0E-10,sdlims(ik,ir,1)))
c              obj(nobj)%quantity(2) = sdlims(ik,ir,1)
              IF ((ielement.NE.0.AND.
     .             opt%obj_option(ielement).EQ.3)) THEN
c     .            opt%ob_stdgrd.EQ.3) THEN
c              IF (opt%ob_stdgrd.EQ.3) THEN
                IF (MOD(isector,2).EQ.1) THEN
                ELSE
c                 IF (ik.LE.3) obj(nobj)%quantity(1:2) = 0.0
                 IF (zs(ik,ir).LT.-1.75) obj(nobj)%quantity(1:2) = 0.0
                ENDIF
              ENDIF
            ENDIF

            obj(nobj)%nmap        = 0

c...        ...
            id = korpg(ik,ir)
 
            DO i1 = 1, 2
              p1(1,i1)   = DBLE(rvertp(i1,id))*DCOS(DBLE(-0.5*dangle))
              p1(2,i1)   = DBLE(zvertp(i1,id))
              p1(3,i1)   = DBLE(rvertp(i1,id))*DSIN(DBLE(-0.5*dangle))
              p1(1,i1+4) = DBLE(rvertp(i1,id))*DCOS(DBLE(+0.5*dangle))
              p1(2,i1+4) = DBLE(zvertp(i1,id))
              p1(3,i1+4) = DBLE(rvertp(i1,id))*DSIN(DBLE(+0.5*dangle))
c              p1(1,i1)   = DBLE(rvertp(i1,id))
c              p1(2,i1)   = DBLE(zvertp(i1,id))
c              p1(3,i1)   = DBLE(rvertp(i1,id))*DTAN(DBLE(-0.5*dangle))
c              p1(1,i1+4) = DBLE(rvertp(i1,id))
c              p1(2,i1+4) = DBLE(zvertp(i1,id))
c              p1(3,i1+4) = DBLE(rvertp(i1,id))*DTAN(DBLE(+0.5*dangle))
            ENDDO
            DO i1 = 3, 4
              p1(1,i1)   = DBLE(rvertp(i1,id))*DCOS(DBLE(-0.5*dangle))
              p1(2,i1)   = DBLE(zvertp(i1,id))
              p1(3,i1)   = DBLE(rvertp(i1,id))*DSIN(DBLE(-0.5*dangle))
              p1(1,i1+4) = DBLE(rvertp(i1,id))*DCOS(DBLE(+0.5*dangle))
              p1(2,i1+4) = DBLE(zvertp(i1,id))
              p1(3,i1+4) = DBLE(rvertp(i1,id))*DSIN(DBLE(+0.5*dangle))
c              p1(1,i1)   = DBLE(rvertp(i1,id))
c              p1(2,i1)   = DBLE(zvertp(i1,id))
c              p1(3,i1)   = DBLE(rvertp(i1,id))*DTAN(DBLE(-0.5*dangle))
c              p1(1,i1+4) = DBLE(rvertp(i1,id))
c              p1(2,i1+4) = DBLE(zvertp(i1,id))
c              p1(3,i1+4) = DBLE(rvertp(i1,id))*DTAN(DBLE(+0.5*dangle))
            ENDDO

            DO i1 = 1, obj(nobj)%nver
              obj(nobj)%v(1:3,i1) = p1(1:3,i1)
            ENDDO
c...        Rotate vertices:
            DO i1 = 1, obj(nobj)%nver
              x1 = obj(nobj)%v(1,i1)
              z1 = obj(nobj)%v(3,i1)
              obj(nobj)%v(1,i1) = DCOS(DBLE(ang)) * x1 -
     .                            DSIN(DBLE(ang)) * z1
              obj(nobj)%v(3,i1) = DSIN(DBLE(ang)) * x1 +
     .                            DCOS(DBLE(ang)) * z1
            ENDDO

            IF ((ielement.NE.0.AND.
     .           (opt%obj_option(ielement).EQ.2.OR.
     .            opt%obj_option(ielement).EQ.4))) THEN
c     .          opt%ob_stdgrd.EQ.2) THEN

c...          Toroidally continuous:

              IF     (ik.EQ.1.AND.ir.GE.irsep) THEN     ! *** THIS IS NOT TRUE IF TARGET TRANSPARENT!
                IF (ielement.NE.0.AND.
     .              opt%obj_option(ielement).EQ.4) THEN
                  obj(nobj)%tsur(1) = SP_GRID_BOUNDARY
                  obj(nobj)%reflec(1) = 0
                ELSE
                  obj(nobj)%tsur(1) = SP_VESSEL_WALL  
                  IF (ielement.NE.0) THEN
                    obj(nobj)%reflec(1) = opt%obj_reflec(ielement)
                  ELSE
c                    obj(nobj)%reflec(1) = opt%ob_stdgrd_reflec
                  ENDIF
                ENDIF
                obj(nobj)%nmap(1) = 1
                obj(nobj)%imap(1,1) = nobj
                obj(nobj)%isur(1,1) = 1
              ELSEIF (ik.EQ.1.AND.ir.LT.irsep) THEN
                obj(nobj)%tsur(1) = SP_GRID_SURFACE
                obj(nobj)%reflec(1) = 0
                obj(nobj)%nmap(1) = 1
                obj(nobj)%imap(1,1) = nobj + cntcell - 1
                obj(nobj)%isur(1,1) = 3
              ELSE
                obj(nobj)%tsur(1) = SP_GRID_SURFACE
                obj(nobj)%reflec(1) = 0
                obj(nobj)%nmap(1) = 1
                obj(nobj)%imap(1,1) = nobj - 1
                obj(nobj)%isur(1,1) = 3
              ENDIF
              obj(nobj)%npts(1) = 2
              obj(nobj)%ipts(1,1) = 2  
              obj(nobj)%ipts(2,1) = 1

              IF (idring(irouts(ik,ir)).EQ.BOUNDARY) THEN
                obj(nobj)%tsur(2) = SP_GRID_BOUNDARY
                obj(nobj)%reflec(2) = 0
c                IF (ielement.NE.0) THEN
c                  obj(nobj)%reflec(2) = opt%obj_reflec(ielement)
c                ELSE
c                  obj(nobj)%reflec(2) = opt%ob_stdgrd_reflec
c                ENDIF
                obj(nobj)%nmap(2) = 1
                obj(nobj)%imap(1,2) = nobj
                obj(nobj)%isur(1,2) = 2
              ELSE
                obj(nobj)%tsur(2) = SP_GRID_SURFACE
                obj(nobj)%reflec(2) = 0
                obj(nobj)%nmap(2) = -MIN(20,noutmap(ik,ir))
                DO i1 = 1, MIN(20,noutmap(ik,ir))
                  obj(nobj)%imap(i1,2) = ikoutmap(ik,ir,i1) 
                  obj(nobj)%isur(i1,2) = iroutmap(ik,ir,i1)
                ENDDO
              ENDIF
              obj(nobj)%npts(2) = 2
              obj(nobj)%ipts(1,2) = 3
              obj(nobj)%ipts(2,2) = 2

              IF     (ik.EQ.ike.AND.ir.GE.irsep) THEN   
                IF (ielement.NE.0.AND.
     .              opt%obj_option(ielement).EQ.4) THEN
                  obj(nobj)%tsur(3) = SP_GRID_BOUNDARY
                  obj(nobj)%reflec(3) = 0 
                ELSE
                  obj(nobj)%tsur(3) = SP_VESSEL_WALL  
                  IF (ielement.NE.0) THEN
                    obj(nobj)%reflec(3) = opt%obj_reflec(ielement)
                  ELSE
c                    obj(nobj)%reflec(3) = opt%ob_stdgrd_reflec
                  ENDIF
                ENDIF
                obj(nobj)%nmap(3) = 1
                obj(nobj)%imap(1,3) = nobj
                obj(nobj)%isur(1,3) = 3
              ELSEIF (ik.EQ.ike.AND.ir.LT.irsep) THEN
                obj(nobj)%tsur(3) = SP_GRID_SURFACE
                obj(nobj)%reflec(3) = 0
                obj(nobj)%nmap(3) = 1
                obj(nobj)%imap(1,3) = nobj - cntcell + 1
                obj(nobj)%isur(1,3) = 1
              ELSE
                obj(nobj)%tsur(3) = SP_GRID_SURFACE
                obj(nobj)%reflec(3) = 0
                obj(nobj)%nmap(3) = 1
                obj(nobj)%imap(1,3) = nobj + 1
                obj(nobj)%isur(1,3) = 1
              ENDIF
              obj(nobj)%npts(3) = 2
              obj(nobj)%ipts(1,3) = 4   
              obj(nobj)%ipts(2,3) = 3

              IF (idring(irins(ik,ir)).EQ.BOUNDARY) THEN
                obj(nobj)%tsur(4) = SP_GRID_BOUNDARY
                obj(nobj)%reflec(4) = 0
                obj(nobj)%nmap(4) = 1
                obj(nobj)%imap(1,4) = nobj
                obj(nobj)%isur(1,4) = 4
              ELSE
                obj(nobj)%tsur(4) = SP_GRID_SURFACE
                obj(nobj)%reflec(4) = 0
                obj(nobj)%nmap(4) = -MIN(20,ninmap(ik,ir))
                DO i1 = 1, MIN(20,ninmap(ik,ir))
                  obj(nobj)%imap(i1,4) = ikinmap(ik,ir,i1) 
                  obj(nobj)%isur(i1,4) = irinmap(ik,ir,i1)
                ENDDO
              ENDIF
              obj(nobj)%npts(4) = 2
              obj(nobj)%ipts(1,4) = 1
              obj(nobj)%ipts(2,4) = 4


            ELSE
c...          Toroidally descretized:
c...
              IF (isector.EQ.1) THEN   
                IF (nsector.EQ.opt%obj_nsector) THEN
                  obj(nobj)%tsur(1) = SP_GRID_SURFACE
                  obj(nobj)%reflec(1) = 0
                  obj(nobj)%nmap(1) = 1
                  obj(nobj)%imap(1,1) = nobj + ncells * (nsector - 1) 
                  obj(nobj)%isur(1,1) = 6 
                ELSE
                  obj(nobj)%tsur(1) = SP_GRID_BOUNDARY  
                  obj(nobj)%reflec(1) = 0
                  obj(nobj)%nmap(1) = 1
                  obj(nobj)%imap(1,1) = nobj 
                  obj(nobj)%isur(1,1) = 1   
                ENDIF
              ELSE
                obj(nobj)%tsur(1) = SP_GRID_SURFACE
                obj(nobj)%reflec(1) = 0
                obj(nobj)%nmap(1) = 1
                obj(nobj)%imap(1,1) = nobj - ncells
                obj(nobj)%isur(1,1) = 6
              ENDIF
              obj(nobj)%npts(1) = 4
              obj(nobj)%ipts(1,1) = 1
              obj(nobj)%ipts(2,1) = 2
              obj(nobj)%ipts(3,1) = 3
              obj(nobj)%ipts(4,1) = 4

              IF     (ik.EQ.1.AND.ir.GE.irsep) THEN     ! *** THIS IS NOT TRUE IF TARGET TRANSPARENT!
                IF (ielement.NE.0.AND.
     .              opt%obj_option(ielement).EQ.5) THEN
                  obj(nobj)%tsur(2) = SP_GRID_BOUNDARY
                  obj(nobj)%reflec(2) = 0 
                ELSE
                  obj(nobj)%tsur(2) = SP_VESSEL_WALL  
                  IF (ielement.NE.0) THEN
                    obj(nobj)%reflec(2) = opt%obj_reflec(ielement)
                  ELSE
c                    obj(nobj)%reflec(2) = opt%ob_stdgrd_reflec
                  ENDIF
                ENDIF
                obj(nobj)%nmap(2) = 1
                obj(nobj)%imap(1,2) = nobj
                obj(nobj)%isur(1,2) = 2
              ELSEIF (ik.EQ.1.AND.ir.LT.irsep) THEN
                obj(nobj)%tsur(2) = SP_GRID_SURFACE
                obj(nobj)%reflec(2) = 0
                obj(nobj)%nmap(2) = 1
                obj(nobj)%imap(1,2) = nobj + cntcell - 1
c                obj(nobj)%imap(1,2) = nobj + ike - 1
                obj(nobj)%isur(1,2) = 4
              ELSE
                obj(nobj)%tsur(2) = SP_GRID_SURFACE
                obj(nobj)%reflec(2) = 0
                obj(nobj)%nmap(2) = 1
                obj(nobj)%imap(1,2) = nobj - 1
                obj(nobj)%isur(1,2) = 4
              ENDIF
              obj(nobj)%npts(2) = 4
c              obj(nobj)%ipts(1,2) = 6
c              obj(nobj)%ipts(2,2) = 5
c              obj(nobj)%ipts(3,2) = 1
c              obj(nobj)%ipts(4,2) = 2
              obj(nobj)%ipts(1,2) = 2   ! If this is right, then the rotation is backwards! 
              obj(nobj)%ipts(2,2) = 1
              obj(nobj)%ipts(3,2) = 5
              obj(nobj)%ipts(4,2) = 6


              IF (idring(irouts(ik,ir)).EQ.BOUNDARY) THEN
                obj(nobj)%tsur(3) = SP_GRID_BOUNDARY
                obj(nobj)%reflec(3) = 0
                obj(nobj)%nmap(3) = 1
                obj(nobj)%imap(1,3) = nobj
                obj(nobj)%isur(1,3) = 3
              ELSE
                obj(nobj)%tsur(3) = SP_GRID_SURFACE
                obj(nobj)%reflec(3) = 0
                obj(nobj)%nmap(3) = -MIN(20,noutmap(ik,ir))
                DO i1 = 1, MIN(20,noutmap(ik,ir))
                  obj(nobj)%imap(i1,3) = ikoutmap(ik,ir,i1) 
                  obj(nobj)%isur(i1,3) = iroutmap(ik,ir,i1)
                ENDDO
              ENDIF
              obj(nobj)%npts(3) = 4
              obj(nobj)%ipts(1,3) = 3
              obj(nobj)%ipts(2,3) = 2
              obj(nobj)%ipts(3,3) = 6
              obj(nobj)%ipts(4,3) = 7

              IF     (ik.EQ.ike.AND.ir.GE.irsep) THEN   
                IF (ielement.NE.0.AND.
     .              opt%obj_option(ielement).EQ.5) THEN
                  obj(nobj)%tsur(4) = SP_GRID_BOUNDARY
                  obj(nobj)%reflec(4) = 0 
                ELSE
                  obj(nobj)%tsur(4) = SP_VESSEL_WALL  
                  IF (ielement.NE.0) THEN
                    obj(nobj)%reflec(4) = opt%obj_reflec(ielement)
                  ELSE
c                    obj(nobj)%reflec(4) = opt%ob_stdgrd_reflec
                  ENDIF
                ENDIF
                obj(nobj)%nmap(4) = 1
                obj(nobj)%imap(1,4) = nobj
                obj(nobj)%isur(1,4) = 4
              ELSEIF (ik.EQ.ike.AND.ir.LT.irsep) THEN
                obj(nobj)%tsur(4) = SP_GRID_SURFACE
                obj(nobj)%reflec(4) = 0
                obj(nobj)%nmap(4) = 1
                obj(nobj)%imap(1,4) = nobj - cntcell + 1
c                obj(nobj)%imap(1,4) = nobj - ike + 1
                obj(nobj)%isur(1,4) = 2
              ELSE
                obj(nobj)%tsur(4) = SP_GRID_SURFACE
                obj(nobj)%reflec(4) = 0
                obj(nobj)%nmap(4) = 1
                obj(nobj)%imap(1,4) = nobj + 1
                obj(nobj)%isur(1,4) = 2
              ENDIF
              obj(nobj)%npts(4) = 4
c              obj(nobj)%ipts(1,4) = 8
c              obj(nobj)%ipts(2,4) = 7
c              obj(nobj)%ipts(3,4) = 3
c              obj(nobj)%ipts(4,4) = 4
              obj(nobj)%ipts(1,4) = 4    ! If this is right, then the rotation is backwards! 
              obj(nobj)%ipts(2,4) = 3
              obj(nobj)%ipts(3,4) = 7
              obj(nobj)%ipts(4,4) = 8

              IF (idring(irins(ik,ir)).EQ.BOUNDARY) THEN
                obj(nobj)%tsur(5) = SP_GRID_BOUNDARY
                obj(nobj)%reflec(5) = 0
                obj(nobj)%nmap(5) = 1
                obj(nobj)%imap(1,5) = nobj
                obj(nobj)%isur(1,5) = 5
              ELSE
                obj(nobj)%tsur(5) = SP_GRID_SURFACE
                obj(nobj)%reflec(5) = 0
                obj(nobj)%nmap(5) = -MIN(20,ninmap(ik,ir))
                DO i1 = 1, MIN(20,ninmap(ik,ir))
                  obj(nobj)%imap(i1,5) = ikinmap(ik,ir,i1) 
                  obj(nobj)%isur(i1,5) = irinmap(ik,ir,i1)
                ENDDO
              ENDIF
              obj(nobj)%npts(5) = 4
              obj(nobj)%ipts(1,5) = 1
              obj(nobj)%ipts(2,5) = 4
              obj(nobj)%ipts(3,5) = 8
              obj(nobj)%ipts(4,5) = 5

              IF (isector.EQ.nsector) THEN 
                IF (nsector.EQ.opt%obj_nsector) THEN
                  obj(nobj)%tsur(6) = SP_GRID_SURFACE
                  obj(nobj)%reflec(6) = 0
                  obj(nobj)%nmap(6) = 1
                  obj(nobj)%imap(1,6) = nobj - ncells * (nsector - 1) 
                  obj(nobj)%isur(1,6) = 1
                ELSE
                  obj(nobj)%tsur(6) = SP_GRID_BOUNDARY
                  obj(nobj)%reflec(6) = 0
                  obj(nobj)%nmap(6) = 1
                  obj(nobj)%imap(1,6) = nobj  
                  obj(nobj)%isur(1,6) = 6  
                ENDIF
              ELSE
                obj(nobj)%tsur(6) = SP_GRID_SURFACE
                obj(nobj)%reflec(6) = 0
                obj(nobj)%nmap(6) = 1
                obj(nobj)%imap(1,6) = nobj + ncells
                obj(nobj)%isur(1,6) = 1
              ENDIF
              obj(nobj)%npts(6) = 4
              obj(nobj)%ipts(1,6) = 8
              obj(nobj)%ipts(2,6) = 7
              obj(nobj)%ipts(3,6) = 6
              obj(nobj)%ipts(4,6) = 5

            ENDIF ! End of GSUR IF block

          ENDDO ! IK loop
        ENDDO ! IR loop

        indmap = 0
        DO i1 = nobj - ncells + 1, nobj
         indmap(obj(i1)%ik,obj(i1)%ir) = i1
        ENDDO

c...    Assign radial cross-field mapping:
        DO i1 = nobj - ncells + 1, nobj
          DO i2 = 2, 4, 2              
            i4 = i2
            IF (nsector.GT.1) i4 = i4 + 1   ! %GSUR=TD_TC if NSECTOR=1 
            IF (obj(i1)%nmap(i4).LT.0) THEN
              DO i3 = 1, -obj(i1)%nmap(i4)
                ik = obj(i1)%imap(i3,i4)
                ir = obj(i1)%isur(i3,i4)
                obj(i1)%imap(i3,i4) = indmap(ik,ir)
                IF     (i4.EQ.2) THEN       ! Almost incomprehensible, even for me, but
                  obj(i1)%isur(i3,i4) = 4   ! these first 2 are for toroidally continuous
                ELSEIF (i4.EQ.4) THEN       ! objects, and the second 2 aren't...
                  obj(i1)%isur(i3,i4) = 2
                ELSEIF (i4.EQ.3) THEN
                  obj(i1)%isur(i3,i4) = 5
                ELSEIF (i4.EQ.5) THEN
                  obj(i1)%isur(i3,i4) = 3
                ELSE
                  CALL ER('ProcessMagneticGrid','Something bad',*99)
                ENDIF
              ENDDO
              obj(i1)%nmap(i4) = -obj(i1)%nmap(i4)
            ENDIF 
          ENDDO
        ENDDO

      ENDDO ! ISECTOR loop

      DEALLOCATE(ninmap)
      DEALLOCATE(ikinmap)
      DEALLOCATE(irinmap)
      DEALLOCATE(noutmap)
      DEALLOCATE(ikoutmap)
      DEALLOCATE(iroutmap)
      DEALLOCATE(indmap)

 98   RETURN
 99   WRITE(0,*) 'IK,IR,IR1 =',ik,ir,ir1
      WRITE(0,*) 'IDRING(IR)=',idring(ir),BOUNDARY
      STOP
      END
c
c ====================================================================== 
c
c
c             side 3
c         3------------4            1  2  3  4  5 ....   nxbin
c         |            |            nxbin+1 ......     2*nxbin
c         |            |            .... 
c         |            |            ....           nybin*nxbin
c  side 2 |            | side 4
c         |            |
c         |            |
c         |            |
c         2------------1
c             side 1
c
c
c 
c
c
      SUBROUTINE LoadInversion(ielement)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

      INTEGER ielement

      INTEGER i1,i2,i3,ix,iy,ncells,isector,nsector,nxbin,nybin,ivol,
     .        istart,iend,fobj,iobj,isid,iobj2,ik,ir
      LOGICAL outofgrid
      REAL    angle,dangle,ang
      REAL*8  p1(3,8),p2(3,8),xcen,ycen,xwidth,ywidth,xdelta,ydelta,
     .        x1,z1,xorigin,yorigin,distsep,distwal,
     .        frac,maxdist,maxdiststd,maxdistxpt


c...    Number of objects per toroidal segment:

      WRITE(0,*) 'LOADING INVERSION'


      nsector = opt%obj_nsector



      fobj = nobj + 1

      IF (ielement.NE.0) THEN
        nxbin = opt%obj_n(ielement,1)
        nybin = opt%obj_n(ielement,2)
        ncells = nxbin * nybin

        xcen = 0.5D0 * (opt%obj_r(ielement,1) + opt%obj_r(ielement,2))
        ycen = 0.5D0 * (opt%obj_z(ielement,1) + opt%obj_z(ielement,2))

        xwidth = opt%obj_r(ielement,2) - opt%obj_r(ielement,1)
        ywidth = opt%obj_z(ielement,2) - opt%obj_z(ielement,1)

        xdelta = xwidth / DBLE(nxbin)
        ydelta = ywidth / DBLE(nybin) 

        xorigin = xcen - 0.5D0 * DBLE(nxbin) * xdelta
        yorigin = ycen - 0.5D0 * DBLE(nybin) * ydelta
      ELSE
      ENDIF


      IF (.TRUE.) THEN
c...    (Perfect) toroidal symmetry:

        ivol = 0  ! Integration volume index 


        maxdiststd = 0.15D0
        maxdistxpt = 0.25D0

        DO iy = nybin, 1, -1
          DO ix = 1, nxbin

            IF (opt%obj_option(ielement).EQ.4) THEN

              xcen = xorigin + xdelta * (DBLE(ix) - 0.5D0)
              ycen = yorigin + ydelta * (DBLE(iy) - 0.5D0)           

              IF (.TRUE.) THEN

              ELSEIF (.FALSE.) THEN
                ik = 1

              ELSE


              ENDIF

            ENDIF


            IF (nobj+1.GT.MAX3D) 
     .        CALL ER('BuildObjects','Insufficient array bounds '//
     .                'for all objects F',*98)     

            ivol = ivol + 1

            nobj = nobj + 1
            obj(nobj)%index       = ielement  ! nobj
            obj(nobj)%type        = OP_INTEGRATION_VOLUME
            obj(nobj)%subtype     = OP_INVERSION_GRID
            obj(nobj)%mode        = 0      
            obj(nobj)%surface     = 1      ! SOLID
            obj(nobj)%wedge1      = 0
            obj(nobj)%wedge2      = 0
            obj(nobj)%colour      = 3
            obj(nobj)%orientation = 1      ! CW
            obj(nobj)%ik          = ix
            obj(nobj)%ir          = iy
            obj(nobj)%in          = 0
            obj(nobj)%ivolume     = ivol
            obj(nobj)%nsur        = 4
            obj(nobj)%gsur(1:4)   = GT_TC
            obj(nobj)%nver        = 4
            obj(nobj)%reflec(1:4) = 0
 


            obj(nobj)%quantity = 1.0
c            obj(nobj)%quantity(1) = 1.0

c            IF (MOD(ix,2).EQ.0) THEN
c              IF (MOD(iy,2).EQ.0) THEN 
c                obj(nobj)%quantity(1) = 1.0
c              ELSE
c                obj(nobj)%quantity(1) = 0.0
c              ENDIF
c            ELSE
c              IF (MOD(iy,2).EQ.0) THEN 
c                obj(nobj)%quantity(1) = 0.0
c              ELSE
c                obj(nobj)%quantity(1) = 1.0
c              ENDIF
c            ENDIF

c            IF (ix.GT.nxbin/2) obj(nobj)%quantity(1) = 0.0

c            obj(3)%quantity(1) = 10.0
c            obj(nobj)%quantity(1) = REAL(ivol)

c...        Vertices:
            obj(nobj)%v(1,1) = xorigin + xdelta * DBLE(ix)   
            obj(nobj)%v(2,1) = yorigin + ydelta * DBLE(iy-1)
            obj(nobj)%v(3,1) = 0.0D0
            obj(nobj)%v(1,2) = xorigin + xdelta * DBLE(ix-1)
            obj(nobj)%v(2,2) = yorigin + ydelta * DBLE(iy-1)
            obj(nobj)%v(3,2) = 0.0D0
            obj(nobj)%v(1,3) = xorigin + xdelta * DBLE(ix-1)
            obj(nobj)%v(2,3) = yorigin + ydelta * DBLE(iy)
            obj(nobj)%v(3,3) = 0.0D0
            obj(nobj)%v(1,4) = xorigin + xdelta * DBLE(ix)
            obj(nobj)%v(2,4) = yorigin + ydelta * DBLE(iy)
            obj(nobj)%v(3,4) = 0.0D0



c...        Surface 1:
            IF (iy.EQ.1) THEN    
              obj(nobj)%tsur(1) = SP_GRID_BOUNDARY
              obj(nobj)%nmap(1) = 1
              obj(nobj)%imap(1,1) = nobj
              obj(nobj)%isur(1,1) = 1
            ELSE
              obj(nobj)%tsur(1) = SP_GRID_SURFACE
              obj(nobj)%nmap(1) = 1
              obj(nobj)%imap(1,1) = nobj + nxbin
              obj(nobj)%isur(1,1) = 3
            ENDIF

            obj(nobj)%npts(1) = 2
            obj(nobj)%ipts(1,1) = 1
            obj(nobj)%ipts(2,1) = 2
c...        Surface 2:
            IF (ix.EQ.1) THEN
              obj(nobj)%tsur(2) = SP_GRID_BOUNDARY
              obj(nobj)%nmap(2) = 1
              obj(nobj)%imap(1,2) = nobj
              obj(nobj)%isur(1,2) = 2
            ELSE
              obj(nobj)%tsur(2) = SP_GRID_SURFACE
              obj(nobj)%nmap(2) = 1
              obj(nobj)%imap(1,2) = nobj - 1
              obj(nobj)%isur(1,2) = 4
            ENDIF
            obj(nobj)%npts(2) = 2
            obj(nobj)%ipts(1,2) = 2
            obj(nobj)%ipts(2,2) = 3
c...        Surface 4:
            IF (iy.EQ.nybin) THEN   
              obj(nobj)%tsur(3) = SP_GRID_BOUNDARY
              obj(nobj)%nmap(3) = 1
              obj(nobj)%imap(1,3) = nobj
              obj(nobj)%isur(1,3) = 3
            ELSE
              obj(nobj)%tsur(3) = SP_GRID_SURFACE
              obj(nobj)%nmap(3) = 1
              obj(nobj)%imap(1,3) = nobj - nxbin
              obj(nobj)%isur(1,3) = 1
            ENDIF
            obj(nobj)%npts(3) = 2
            obj(nobj)%ipts(1,3) = 3
            obj(nobj)%ipts(2,3) = 4
c...        Surface 4:
            IF (ix.EQ.nxbin) THEN
              obj(nobj)%tsur(4) = SP_GRID_BOUNDARY
              obj(nobj)%nmap(4) = 1
              obj(nobj)%imap(1,4) = nobj
              obj(nobj)%isur(1,4) = 4
            ELSE
              obj(nobj)%tsur(4) = SP_GRID_SURFACE
              obj(nobj)%nmap(4) = 1
              obj(nobj)%imap(1,4) = nobj + 1
              obj(nobj)%isur(1,4) = 2
            ENDIF
            obj(nobj)%npts(4) = 2
            obj(nobj)%ipts(1,4) = 4
            obj(nobj)%ipts(2,4) = 1

          ENDDO
        ENDDO



        IF (opt%obj_option(ielement).EQ.4) THEN
c...      Need to build connection map since grid is no longer regular:

          DO iobj = fobj, nobj
            DO isid = 1, obj(iobj)%nsur
 
              SELECTCASE (isid)
                CASE(1)
                  obj(iobj)%tsur(isid) = SP_GRID_BOUNDARY
                  obj(iobj)%imap(1,isid) = iobj
                  obj(iobj)%isur(1,isid) = 1                                     
                  DO iobj2 = MIN(nobj,iobj+1), MIN(nobj,iobj+nxbin+1) 
                    IF (obj(iobj)%v(1,1).EQ.obj(iobj2)%v(1,4).AND.
     .                  obj(iobj)%v(2,1).EQ.obj(iobj2)%v(2,4).AND.
     .                  obj(iobj)%v(1,2).EQ.obj(iobj2)%v(1,3).AND.
     .                  obj(iobj)%v(2,2).EQ.obj(iobj2)%v(2,3)) THEN
                      obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                      obj(iobj)%imap(1,isid) = iobj2
                      obj(iobj)%isur(1,isid) = 3                   
                      EXIT
                    ENDIF
                  ENDDO

                CASE(2)
                  iobj2 = MAX(fobj,iobj-1)
                  IF (obj(iobj)%v(1,2).EQ.obj(iobj2)%v(1,1).AND.
     .                obj(iobj)%v(2,2).EQ.obj(iobj2)%v(2,1).AND.
     .                obj(iobj)%v(1,3).EQ.obj(iobj2)%v(1,4).AND.
     .                obj(iobj)%v(2,3).EQ.obj(iobj2)%v(2,4)) THEN
                    obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                    obj(iobj)%imap(1,isid) = iobj2
                    obj(iobj)%isur(1,isid) = 4                    
                  ELSE
                    obj(iobj)%tsur(isid) = SP_GRID_BOUNDARY
                    obj(iobj)%imap(1,isid) = iobj
                    obj(iobj)%isur(1,isid) = 2                                      
                  ENDIF

                CASE(3)
                  obj(iobj)%tsur(isid) = SP_GRID_BOUNDARY
                  obj(iobj)%imap(1,isid) = iobj
                  obj(iobj)%isur(1,isid) = 3                                    
                  DO iobj2 = MAX(fobj,iobj-1),MAX(fobj,iobj-nxbin-1),-1
                    IF (obj(iobj)%v(1,3).EQ.obj(iobj2)%v(1,2).AND.
     .                  obj(iobj)%v(2,3).EQ.obj(iobj2)%v(2,2).AND.
     .                  obj(iobj)%v(1,4).EQ.obj(iobj2)%v(1,1).AND.
     .                  obj(iobj)%v(2,4).EQ.obj(iobj2)%v(2,1)) THEN
                      obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                      obj(iobj)%imap(1,isid) = iobj2
                      obj(iobj)%isur(1,isid) = 1                  
                      EXIT
                    ENDIF
                  ENDDO               

                CASE(4)
                  iobj2 = MIN(nobj,iobj+1)
                  IF (obj(iobj)%v(1,1).EQ.obj(iobj2)%v(1,2).AND.
     .                obj(iobj)%v(2,1).EQ.obj(iobj2)%v(2,2).AND.
     .                obj(iobj)%v(1,4).EQ.obj(iobj2)%v(1,3).AND.
     .                obj(iobj)%v(2,4).EQ.obj(iobj2)%v(2,3)) THEN
                    obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                    obj(iobj)%imap(1,isid) = iobj2
                    obj(iobj)%isur(1,isid) = 2                   
                  ELSE
                    obj(iobj)%tsur(isid) = SP_GRID_BOUNDARY
                    obj(iobj)%imap(1,isid) = iobj
                    obj(iobj)%isur(1,isid) = 4                                    
                  ENDIF

                CASEDEFAULT
              ENDSELECT       

            ENDDO
 
c            WRITE(0,*) 'CEL:',iobj,obj(iobj)%tsur(1:4)
c            WRITE(0,*) '   :',iobj,obj(iobj)%imap(1,1:4)
c            WRITE(0,*) '   :',iobj,obj(iobj)%isur(1,1:4)

          ENDDO

        ENDIF


      ELSE
        CALL ER('LoadInversion','Unknown option',*99)
      ENDIF

c...  Local mesh refinement adjustment based on mask, from a previous 
c     iteration:


c...  Generalized cell mapping method:

      WRITE(0,*) 'DONE'


 98   RETURN
 99   STOP
      END


c
c ======================================================================
c
c subroutine: LoadVesselStructures
c
c
      SUBROUTINE LoadVesselStructures(ielement)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

c...  Input:
      INTEGER ielement

      INTEGER AddVertex,AddSurface

      INTEGER fp,count,i1,idum1,istart,iobj,opt_order,opt_flip(3)
      REAL*8  newvtx(3,3),angle,mat(3,3)               ! Assumes triangles from .raw file!
      CHARACTER file*1024,buffer*2048
      TYPE(type_surface) newsrf


c...unit scaling, rotation, toroidal sub-section, etc.
c      opt%ob_raw_num = 2
c      opt%ob_raw_ind(1,1) = 1
c      opt%ob_raw_ind(1,2) = 100
c      opt%ob_raw_colour(1) = 1
c      opt%ob_raw_fname(1) = 'P3-upper.raw'

c      opt%ob_raw_ind(2,1) = 1
c      opt%ob_raw_ind(2,2) = 10000
c      opt%ob_raw_colour(2) = 1
c      opt%ob_raw_fname(2) = 'HU01-rogowski.raw'


c...   opt%obj_orientation = abcd
c          a - flip about x-axis
c          b - flip about y-axis
c          c - flip about z-axis
c          d - change clockwise ordering of points

      opt_order = 1
      opt_flip  = 0
      opt_order    = MOD(opt%obj_orientation(ielement)     ,10)
      opt_flip(1)  = MOD(opt%obj_orientation(ielement)/1000,10)
      opt_flip(2)  = MOD(opt%obj_orientation(ielement)/100 ,10)
      opt_flip(3)  = MOD(opt%obj_orientation(ielement)/10  ,10)
      WRITE(0,*) 'ORDER,FLIP=',opt_order,opt_flip


      fp = 99


       WRITE(0,*) 'LOADING VESSEL STRUCTURES!'//
     .     opt%obj_fname(ielement)

c...  Load up surface and vertex arrays:

      istart = nsrf + 1

c      file = 'HU08-port-plate.raw'
c      file = 'test-p3.raw'
c      file = opt%ob_raw_fname(iraw)
      file = opt%obj_fname(ielement)

c...  Rotate about y-axis (swing):
      IF (opt%obj_yangle(ielement).NE.0.0) THEN
        CALL Calc_Transform2(mat,0.0D0,1,0)
        angle = DBLE(opt%obj_yangle(ielement)*3.141592/180.0)
        CALL Calc_Transform2(mat,angle,2,1)
      ENDIF

      OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .     FORM='FORMATTED',STATUS='OLD',ERR=98)     
 
      count = 0
      DO WHILE (.TRUE.)       

        READ(fp,'(A2048)',END=10) buffer         
        IF (buffer(1:3).EQ.'obj') THEN
          READ(buffer(5:),*) iobj
c          WRITE(0,*) 'OBJNUM:',iobj
          DO WHILE (.TRUE.) 
            READ(fp,'(A2048)',END=10) buffer 
            IF (buffer(1:3).NE.'obj') THEN
              IF ((iobj.GE.opt%obj_n(ielement,1).OR.
     .             -1  .EQ.opt%obj_n(ielement,1)).AND.
     .            (iobj.LE.opt%obj_n(ielement,2).OR.
     .             -1  .EQ.opt%obj_n(ielement,2))) THEN
c...            Extract vertices:

c               DO program to RAY convention: x=x, y=-z, z=y

                SELECTCASE (opt_order)
                  CASE(1)
                    READ(buffer,*) 
     .                (newvtx(1,i1),newvtx(3,i1),newvtx(2,i1),i1=1,3)
                  CASE(2)
                    READ(buffer,*) 
     .                (newvtx(1,i1),newvtx(3,i1),newvtx(2,i1),i1=3,1,-1)
                  CASE DEFAULT
                    CALL ER('LoadVesselStructures','Unknown '//
     .                      'orientation option',*99)
                ENDSELECT

c...            Flipping (in case orientation here not consisitent with source program): 
                DO i1 = 1, 3
                  IF (opt_flip(i1).EQ.1) newvtx(i1,1:3)= -newvtx(i1,1:3)
                ENDDO

                DO i1 = 1, 3
c...             Convert units from mm to m:
                 newvtx(1:3,i1)=newvtx(1:3,i1)*opt%obj_scale(ielement)  !0.001D0
c...             Toroidal rotation (about y axis):
                 IF (opt%obj_yangle(ielement).NE.0.0) THEN
                   CALL Transform_Vect(mat,newvtx(1,i1))
                 ENDIF
                ENDDO

                newsrf%type = SP_PLANAR_POLYGON
                newsrf%nvtx = 3
                DO i1 = 1, 3
                  newsrf%ivtx(i1) = AddVertex(newvtx(1,i1))
                ENDDO
                idum1 = AddSurface(newsrf)
              ENDIF
            ELSE
c              WRITE(0,*) '  DONE, GETTING OUT OF LOOP',count
              BACKSPACE fp
              EXIT
            ENDIF
          ENDDO
          count = count + 1
        ENDIF
        
      ENDDO
 10   CONTINUE
      CLOSE(fp)

      WRITE(0,*) 'VERTICES? SURFACES?',nvtx,nsrf
      WRITE(0,*) '                   ',istart

c...  Assign object(s):    

      WRITE(0,*) 'NOBJ:',nobj,MAX3D


      IF (nobj+1.GT.MAX3D) 
     .  CALL ER('LoadVesselStructures','Insufficient array bounds '//
     .          'for all objects G',*99)    

      IF (istart.GT.nsrf) THEN
        WRITE(0,*) 'LoadVesselStructures: Strange, no objects loaded'
        RETURN
      ENDIF

      nobj = nobj + 1
      WRITE(0,*) 'VESSEL STRUCTURE IOBJ:',nobj

      obj(nobj)%index       = ielement  ! nobj
      obj(nobj)%type        = OP_EMPTY
      obj(nobj)%mode        = 0      
      obj(nobj)%surface     = 1      ! SOLID
      obj(nobj)%wedge1      = 0
      obj(nobj)%wedge2      = 0
      obj(nobj)%colour      = 1
      obj(nobj)%orientation = 1      ! CW
      obj(nobj)%ik          = 0
      obj(nobj)%ir          = 0
      obj(nobj)%in          = -1  ! What should this be?
      obj(nobj)%ivolume     = 0
      obj(nobj)%nside       = 1
      obj(nobj)%iside(1,1)  = istart ! Start index of range of surfaces in surface array, from loading code above
      obj(nobj)%iside(1,2)  = nsrf   ! End index of range of surfaces in surface array
      obj(nobj)%gsur(1)     = GT_TD
      obj(nobj)%tsur(1)     = SP_VESSEL_WALL
      obj(nobj)%reflec(1)   = opt%obj_reflec(ielement)
c..   Defunct:
      obj(nobj)%nsur        = 0
      obj(nobj)%ipts(2,1)   = 0
      obj(nobj)%nmap(1)     = 0

      WRITE(0,*) 'DONE'




      RETURN
 98   WRITE(0,*) 'ERROR LoadVesselStructures: File not found'
      WRITE(0,*) '   '//file(1:LEN_TRIM(file))
 99   STOP
      END

c
c ====================================================================== 
c
c
c             side 3
c         3------------4            1  2  3  4  5 ....   nxbin
c         |            |            nxbin+1 ......     2*nxbin
c         |            |            .... 
c         |            |            ....           nybin*nxbin
c  side 2 |            | side 4
c         |            |
c         |            |
c         |            |
c         2------------1
c             side 1
c
c
c 
c
c
      SUBROUTINE BuildInversionMesh_Old(ielement)
      USE mod_out985
      USE mod_out985_variables
      use mod_params
      use mod_cgeom
      use mod_slcom
      IMPLICIT none

      INTEGER ielement


c     INCLUDE 'params'
c     INCLUDE 'cgeom'    ! rxp,zxp used below
c      INCLUDE 'comtor'
c      INCLUDE 'pindata'
c      include 'dynam2'
c     INCLUDE 'slcom'    ! eirntorseg used


      LOGICAL CheckInversionCell

      INTEGER i1,i2,i3,ix,iy,ncells,isector,nsector,nxbin,nybin,ivol,
     .        istart,iend,fobj,iobj,isid,iobj2,ik,ir
      LOGICAL outofgrid
      REAL    angle,dangle,ang
      REAL*8  p1(3,8),p2(3,8),xcen,ycen,xwidth,ywidth,xdelta,ydelta,
     .        x1,z1,xorigin,yorigin,distsep,distwal,
     .        frac,maxdist,maxdiststd,maxdistxpt


c...    Number of objects per toroidal segment:

c      WRITE(0,*) 'BUILDING RADIAL MAP'


      nsector = opt%obj_nsector
      IF (nsector.EQ.-1) nsector = eirntorseg

      dangle = 360.0 / REAL(nsector) / RADDEG


      fobj = nobj + 1

      IF (ielement.NE.0) THEN
        nxbin = opt%obj_n(ielement,1)
        nybin = opt%obj_n(ielement,2)
        ncells = nxbin * nybin

        xcen = 0.5D0 * (opt%obj_r(ielement,1) + opt%obj_r(ielement,2))
        ycen = 0.5D0 * (opt%obj_z(ielement,1) + opt%obj_z(ielement,2))

        xwidth = opt%obj_r(ielement,2) - opt%obj_r(ielement,1)
        ywidth = opt%obj_z(ielement,2) - opt%obj_z(ielement,1)

        xdelta = xwidth / DBLE(nxbin)
        ydelta = ywidth / DBLE(nybin) 

        xorigin = xcen - 0.5D0 * DBLE(nxbin) * xdelta
        yorigin = ycen - 0.5D0 * DBLE(nybin) * ydelta
      ELSE
c...    Original code:
        STOP 'NO LONGER SUPPORTED: 98798'
c        nxbin = opt%ob_invgrd_nxbin
c        nybin = opt%ob_invgrd_nybin
        ncells = nxbin * nybin

c        xcen = opt%ob_invgrd_xcen
c        ycen = opt%ob_invgrd_ycen 
        IF (xcen.EQ.-99.0) xcen = rxp
        IF (ycen.EQ.-99.0) ycen = zxp

c        xwidth = opt%ob_invgrd_xwidth ! 0.4
c        ywidth = opt%ob_invgrd_ywidth ! 0.4
        xdelta = xwidth / DBLE(nxbin)
        ydelta = ywidth / DBLE(nybin)
 
        xorigin = xcen - 0.5D0 * DBLE(nxbin) * xdelta
        yorigin = ycen - 0.5D0 * DBLE(nybin) * ydelta
      ENDIF


      IF ((ielement.NE.0.AND.
     .     (opt%obj_option(ielement).EQ.2.OR.
     .      opt%obj_option(ielement).EQ.4))) THEN
c     .    opt%ob_invgrd.EQ.2) THEN
c...    (Perfect) toroidal symmetry:

        ivol = 0  ! Integration volume index 


        maxdiststd = 0.15D0
        maxdistxpt = 0.25D0

        DO iy = nybin, 1, -1
          DO ix = 1, nxbin

            IF (opt%obj_option(ielement).EQ.4) THEN

              xcen = xorigin + xdelta * (DBLE(ix) - 0.5D0)
              ycen = yorigin + ydelta * (DBLE(iy) - 0.5D0)           


              IF (.TRUE.) THEN
                IF (.NOT.CheckInversionCell(1,xcen,ycen)) CYCLE

              ELSEIF (.FALSE.) THEN
                ik = 1
                ir = irsep-1
                CALL GridPos(ik,ir,SNGL(xcen),SNGL(ycen),
     .                       .TRUE.,outofgrid)
                IF (outofgrid) CYCLE

              ELSE

                CALL GetSepDist(2,xcen,ycen,distsep)
                distwal = 1000.0D0
c                CALL GetWalDist(1,xcen,ycen,distwal)

c                IF (iy.EQ.1) WRITE(0,*) 
c     .            'MIND:',distsep,MIN(distsep,distwal)

                IF (nrs.EQ.65) THEN
                  frac = MIN(1.0D0,MIN(DABS(ycen-DBLE(zxp)),
     .                                 DABS(ycen-(-1.2D0) ))/0.35D0)
c                  frac = MIN(1.0D0,MIN(DABS(ycen-DBLE(zxp)),
c     .                                 DABS(ycen-(-1.2D0) ),
c     .                                 DABS(xcen-DBLE(rxp)),
c     .                                 DABS(xcen-  0.6D0  ))/0.35D0)
                  maxdist = (1.0 - frac)*maxdistxpt + frac*maxdiststd
c                  IF (ix.EQ.1) WRITE(0,*) 'WHAT!',iy,ycen,frac
                ELSE
                  maxdist = maxdiststd
                ENDIF
              
                IF (MIN(distsep,distwal).GT.maxdist) CYCLE

              ENDIF

            ENDIF


            IF (nobj+1.GT.MAX3D) 
     .        CALL ER('BuildObjects','Insufficient array bounds '//
     .                'for all objects H',*98)     

            ivol = ivol + 1

            nobj = nobj + 1
            obj(nobj)%index       = ielement  ! nobj
            obj(nobj)%type        = OP_INTEGRATION_VOLUME
            obj(nobj)%subtype     = OP_INVERSION_GRID
            obj(nobj)%mode        = 0      
            obj(nobj)%surface     = 1      ! SOLID
            obj(nobj)%wedge1      = 0
            obj(nobj)%wedge2      = 0
            obj(nobj)%colour      = 3
            obj(nobj)%orientation = 1      ! CW
            obj(nobj)%ik          = ix
            obj(nobj)%ir          = iy
            obj(nobj)%in          = 0
            obj(nobj)%ivolume     = ivol
            obj(nobj)%nsur        = 4
            obj(nobj)%gsur(1:4)   = GT_TC
            obj(nobj)%nver        = 4
            obj(nobj)%reflec(1:4) = 0
 


            obj(nobj)%quantity = 1.0
c            obj(nobj)%quantity(1) = 1.0

c...        Checkerboard:
            IF (MOD(ix,2).EQ.0) THEN
              IF (MOD(iy,2).EQ.0) THEN 
                obj(nobj)%quantity(1) = 1.0
              ELSE
                obj(nobj)%quantity(1) = 0.0
              ENDIF
            ELSE
              IF (MOD(iy,2).EQ.0) THEN 
                obj(nobj)%quantity(1) = 0.0
              ELSE
                obj(nobj)%quantity(1) = 1.0
              ENDIF
            ENDIF

c            IF (ix.GT.nxbin/2) obj(nobj)%quantity(1) = 0.0

c            obj(3)%quantity(1) = 10.0
c            obj(nobj)%quantity(1) = REAL(ivol)

c...        Vertices:
            obj(nobj)%v(1,1) = xorigin + xdelta * DBLE(ix)   
            obj(nobj)%v(2,1) = yorigin + ydelta * DBLE(iy-1)
            obj(nobj)%v(3,1) = 0.0D0
            obj(nobj)%v(1,2) = xorigin + xdelta * DBLE(ix-1)
            obj(nobj)%v(2,2) = yorigin + ydelta * DBLE(iy-1)
            obj(nobj)%v(3,2) = 0.0D0
            obj(nobj)%v(1,3) = xorigin + xdelta * DBLE(ix-1)
            obj(nobj)%v(2,3) = yorigin + ydelta * DBLE(iy)
            obj(nobj)%v(3,3) = 0.0D0
            obj(nobj)%v(1,4) = xorigin + xdelta * DBLE(ix)
            obj(nobj)%v(2,4) = yorigin + ydelta * DBLE(iy)
            obj(nobj)%v(3,4) = 0.0D0



c...        Surface 1:
            IF (iy.EQ.1) THEN    
              obj(nobj)%tsur(1) = SP_GRID_BOUNDARY
              obj(nobj)%nmap(1) = 1
              obj(nobj)%imap(1,1) = nobj
              obj(nobj)%isur(1,1) = 1
            ELSE
              obj(nobj)%tsur(1) = SP_GRID_SURFACE
              obj(nobj)%nmap(1) = 1
              obj(nobj)%imap(1,1) = nobj + nxbin
              obj(nobj)%isur(1,1) = 3
            ENDIF

            obj(nobj)%npts(1) = 2
            obj(nobj)%ipts(1,1) = 1
            obj(nobj)%ipts(2,1) = 2
c...        Surface 2:
            IF (ix.EQ.1) THEN
              obj(nobj)%tsur(2) = SP_GRID_BOUNDARY
              obj(nobj)%nmap(2) = 1
              obj(nobj)%imap(1,2) = nobj
              obj(nobj)%isur(1,2) = 2
            ELSE
              obj(nobj)%tsur(2) = SP_GRID_SURFACE
              obj(nobj)%nmap(2) = 1
              obj(nobj)%imap(1,2) = nobj - 1
              obj(nobj)%isur(1,2) = 4
            ENDIF
            obj(nobj)%npts(2) = 2
            obj(nobj)%ipts(1,2) = 2
            obj(nobj)%ipts(2,2) = 3
c...        Surface 4:
            IF (iy.EQ.nybin) THEN   
              obj(nobj)%tsur(3) = SP_GRID_BOUNDARY
              obj(nobj)%nmap(3) = 1
              obj(nobj)%imap(1,3) = nobj
              obj(nobj)%isur(1,3) = 3
            ELSE
              obj(nobj)%tsur(3) = SP_GRID_SURFACE
              obj(nobj)%nmap(3) = 1
              obj(nobj)%imap(1,3) = nobj - nxbin
              obj(nobj)%isur(1,3) = 1
            ENDIF
            obj(nobj)%npts(3) = 2
            obj(nobj)%ipts(1,3) = 3
            obj(nobj)%ipts(2,3) = 4
c...        Surface 4:
            IF (ix.EQ.nxbin) THEN
              obj(nobj)%tsur(4) = SP_GRID_BOUNDARY
              obj(nobj)%nmap(4) = 1
              obj(nobj)%imap(1,4) = nobj
              obj(nobj)%isur(1,4) = 4
            ELSE
              obj(nobj)%tsur(4) = SP_GRID_SURFACE
              obj(nobj)%nmap(4) = 1
              obj(nobj)%imap(1,4) = nobj + 1
              obj(nobj)%isur(1,4) = 2
            ENDIF
            obj(nobj)%npts(4) = 2
            obj(nobj)%ipts(1,4) = 4
            obj(nobj)%ipts(2,4) = 1

          ENDDO
        ENDDO



        IF (opt%obj_option(ielement).EQ.4) THEN
c...      Need to build connection map since grid is no longer regular:

          DO iobj = fobj, nobj
            DO isid = 1, obj(iobj)%nsur
 
              SELECTCASE (isid)
                CASE(1)
                  obj(iobj)%tsur(isid) = SP_GRID_BOUNDARY
                  obj(iobj)%imap(1,isid) = iobj
                  obj(iobj)%isur(1,isid) = 1                                     
                  DO iobj2 = MIN(nobj,iobj+1), MIN(nobj,iobj+nxbin+1) 
                    IF (obj(iobj)%v(1,1).EQ.obj(iobj2)%v(1,4).AND.
     .                  obj(iobj)%v(2,1).EQ.obj(iobj2)%v(2,4).AND.
     .                  obj(iobj)%v(1,2).EQ.obj(iobj2)%v(1,3).AND.
     .                  obj(iobj)%v(2,2).EQ.obj(iobj2)%v(2,3)) THEN
                      obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                      obj(iobj)%imap(1,isid) = iobj2
                      obj(iobj)%isur(1,isid) = 3                   
                      EXIT
                    ENDIF
                  ENDDO

                CASE(2)
                  iobj2 = MAX(fobj,iobj-1)
                  IF (obj(iobj)%v(1,2).EQ.obj(iobj2)%v(1,1).AND.
     .                obj(iobj)%v(2,2).EQ.obj(iobj2)%v(2,1).AND.
     .                obj(iobj)%v(1,3).EQ.obj(iobj2)%v(1,4).AND.
     .                obj(iobj)%v(2,3).EQ.obj(iobj2)%v(2,4)) THEN
                    obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                    obj(iobj)%imap(1,isid) = iobj2
                    obj(iobj)%isur(1,isid) = 4                    
                  ELSE
                    obj(iobj)%tsur(isid) = SP_GRID_BOUNDARY
                    obj(iobj)%imap(1,isid) = iobj
                    obj(iobj)%isur(1,isid) = 2                                      
                  ENDIF

                CASE(3)
                  obj(iobj)%tsur(isid) = SP_GRID_BOUNDARY
                  obj(iobj)%imap(1,isid) = iobj
                  obj(iobj)%isur(1,isid) = 3                                    
                  DO iobj2 = MAX(fobj,iobj-1),MAX(fobj,iobj-nxbin-1),-1
                    IF (obj(iobj)%v(1,3).EQ.obj(iobj2)%v(1,2).AND.
     .                  obj(iobj)%v(2,3).EQ.obj(iobj2)%v(2,2).AND.
     .                  obj(iobj)%v(1,4).EQ.obj(iobj2)%v(1,1).AND.
     .                  obj(iobj)%v(2,4).EQ.obj(iobj2)%v(2,1)) THEN
                      obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                      obj(iobj)%imap(1,isid) = iobj2
                      obj(iobj)%isur(1,isid) = 1                  
                      EXIT
                    ENDIF
                  ENDDO               

                CASE(4)
                  iobj2 = MIN(nobj,iobj+1)
                  IF (obj(iobj)%v(1,1).EQ.obj(iobj2)%v(1,2).AND.
     .                obj(iobj)%v(2,1).EQ.obj(iobj2)%v(2,2).AND.
     .                obj(iobj)%v(1,4).EQ.obj(iobj2)%v(1,3).AND.
     .                obj(iobj)%v(2,4).EQ.obj(iobj2)%v(2,3)) THEN
                    obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                    obj(iobj)%imap(1,isid) = iobj2
                    obj(iobj)%isur(1,isid) = 2                   
                  ELSE
                    obj(iobj)%tsur(isid) = SP_GRID_BOUNDARY
                    obj(iobj)%imap(1,isid) = iobj
                    obj(iobj)%isur(1,isid) = 4                                    
                  ENDIF

                CASEDEFAULT
              ENDSELECT       

            ENDDO
 
c            WRITE(0,*) 'CEL:',iobj,obj(iobj)%tsur(1:4)
c            WRITE(0,*) '   :',iobj,obj(iobj)%imap(1,1:4)
c            WRITE(0,*) '   :',iobj,obj(iobj)%isur(1,1:4)

          ENDDO

        ENDIF


      ELSEIF ((ielement.NE.0.AND.
     .         opt%obj_option(ielement).EQ.1)) THEN
c     .        opt%ob_invgrd.EQ.1) THEN

c...    Toroidal approximation via discretization:

        istart = 1
        iend   = nsector 
c        istart = -4
c        iend   = nsector / 2

        IF (iend-istart+1.GT.nsector) 
     .    CALL ER('BuildInversionMesh','An excess of toroidal '//
     .            'sectors requested',*99)

        DO isector = istart, iend

          ang = REAL(isector - 1) * dangle
  
          ivol = 0  ! Integration volume index 

          DO iy = nybin, 1, -1
            DO ix = 1, nxbin

              IF (nobj+1.GT.MAX3D) 
     .          CALL ER('BuildObjects','Insufficient array bounds '//
     .                  'for all objects I',*98)     

              ivol = ivol + 1

              nobj = nobj + 1
              obj(nobj)%index       = ielement  ! nobj
              obj(nobj)%type        = OP_INTEGRATION_VOLUME
              obj(nobj)%mode        = 0      
              obj(nobj)%surface     = 1      ! SOLID
              obj(nobj)%phi         = ang
              obj(nobj)%wedge1      = 0
              obj(nobj)%wedge2      = 0
              obj(nobj)%colour      = 3
              obj(nobj)%orientation = 1      ! CW
              obj(nobj)%ik          = ix
              obj(nobj)%ir          = iy
              obj(nobj)%in          = 0
              obj(nobj)%ivolume     = ivol
              obj(nobj)%nsur        = 6
              obj(nobj)%gsur(1:6)   = GT_TD
              obj(nobj)%nver        = 8

c             obj(nobj)%quantity = 0.0
c              obj(nobj)%quantity(1) = 1.0


c              IF (MOD(ix,2).EQ.0) THEN
c                IF (MOD(iy,2).EQ.0) THEN 
c                  obj(nobj)%quantity(1) = 1.0
c                ELSE
c                  obj(nobj)%quantity(1) = 0.0
c                ENDIF
c              ELSE
c                IF (MOD(iy,2).EQ.0) THEN 
c                  obj(nobj)%quantity(1) = 0.0
c                ELSE
c                  obj(nobj)%quantity(1) = 1.0
c                ENDIF
c              ENDIF
c              obj(3)%quantity(1) = 10.0
c            obj(nobj)%quantity(1) = REAL(ivol)
c            obj(nobj)%quantity(2) = 1.0

              obj(nobj)%nmap = 0

c...          ...
              p1(1,1)   = xorigin + xdelta * DBLE(ix)   ! Clean this up...
              p1(2,1)   = yorigin + ydelta * DBLE(iy-1)
              p1(3,1)   = DBLE(p1(1,1))*DTAN(DBLE(-0.5*dangle))
              p1(1,1+4) = p1(1,1)
              p1(2,1+4) = p1(2,1)
              p1(3,1+4) = DBLE(p1(1,1))*DTAN(DBLE(+0.5*dangle))

              p1(1,2)   = xorigin + xdelta * DBLE(ix-1)
              p1(2,2)   = yorigin + ydelta * DBLE(iy-1)
              p1(3,2)   = DBLE(p1(1,2))*DTAN(DBLE(-0.5*dangle))
              p1(1,2+4) = p1(1,2)
              p1(2,2+4) = p1(2,2)
              p1(3,2+4) = DBLE(p1(1,2))*DTAN(DBLE(+0.5*dangle))

              p1(1,3)   = xorigin + xdelta * DBLE(ix-1)
              p1(2,3)   = yorigin + ydelta * DBLE(iy)
              p1(3,3)   = DBLE(p1(1,3))*DTAN(DBLE(-0.5*dangle))
              p1(1,3+4) = p1(1,3)
              p1(2,3+4) = p1(2,3)
              p1(3,3+4) = DBLE(p1(1,3))*DTAN(DBLE(+0.5*dangle))

              p1(1,4)   = xorigin + xdelta * DBLE(ix)
              p1(2,4)   = yorigin + ydelta * DBLE(iy)
              p1(3,4)   = DBLE(p1(1,4))*DTAN(DBLE(-0.5*dangle))
              p1(1,4+4) = p1(1,4)
              p1(2,4+4) = p1(2,4)
              p1(3,4+4) = DBLE(p1(1,4))*DTAN(DBLE(+0.5*dangle))
c... 
              DO i1 = 1, 8
                obj(nobj)%v(1:3,i1) = p1(1:3,i1)
              ENDDO
c...          Rotate vertices:
              DO i1 = 1, 8
                x1 = obj(nobj)%v(1,i1)
                z1 = obj(nobj)%v(3,i1)
                obj(nobj)%v(1,i1) = DCOS(DBLE(ang)) * x1 -
     .                              DSIN(DBLE(ang)) * z1
                obj(nobj)%v(3,i1) = DSIN(DBLE(ang)) * x1 +
     .                              DCOS(DBLE(ang)) * z1
              ENDDO
c...          Surface 1:
              IF (isector.EQ.istart) THEN
                IF (iend-istart+1.EQ.nsector) THEN   
c...              Full torus:
                  obj(nobj)%tsur(1) = SP_GRID_SURFACE
                  obj(nobj)%nmap(1) = 1
                  obj(nobj)%imap(1,1) = nobj + ncells * (nsector - 1)
                  obj(nobj)%isur(1,1) = 6
                ELSE
                  obj(nobj)%tsur(1) = SP_GRID_BOUNDARY
                  obj(nobj)%nmap(1) = 1
                  obj(nobj)%imap(1,1) = nobj
                  obj(nobj)%isur(1,1) = 1
                ENDIF
              ELSE
                obj(nobj)%tsur(1) = SP_GRID_SURFACE
                obj(nobj)%nmap(1) = 1
                obj(nobj)%imap(1,1) = nobj - ncells
                obj(nobj)%isur(1,1) = 6
              ENDIF
              obj(nobj)%npts(1) = 4
              obj(nobj)%ipts(1,1) = 1
              obj(nobj)%ipts(2,1) = 2
              obj(nobj)%ipts(3,1) = 3
              obj(nobj)%ipts(4,1) = 4
c...          Surface 2:
              IF (iy.EQ.1) THEN     ! *** THIS IS NOT TRUE IF TARGET TRANSPARENT!
                obj(nobj)%tsur(2) = SP_GRID_BOUNDARY
                obj(nobj)%nmap(2) = 1
                obj(nobj)%imap(1,2) = nobj
                obj(nobj)%isur(1,2) = 2
              ELSE
                obj(nobj)%tsur(2) = SP_GRID_SURFACE
                obj(nobj)%nmap(2) = 1
                obj(nobj)%imap(1,2) = nobj + nxbin
                obj(nobj)%isur(1,2) = 4
              ENDIF
              obj(nobj)%npts(2) = 4
              obj(nobj)%ipts(1,2) = 2   ! If this is right, then the rotation is backwards! 
              obj(nobj)%ipts(2,2) = 1
              obj(nobj)%ipts(3,2) = 5
              obj(nobj)%ipts(4,2) = 6
c...          Surface 3:
              IF (ix.EQ.1) THEN
                obj(nobj)%tsur(3) = SP_GRID_BOUNDARY
                obj(nobj)%nmap(3) = 1
                obj(nobj)%imap(1,3) = nobj
                obj(nobj)%isur(1,3) = 3
              ELSE
                obj(nobj)%tsur(3) = SP_GRID_SURFACE
                obj(nobj)%nmap(3) = 1
                obj(nobj)%imap(1,3) = nobj - 1
                obj(nobj)%isur(1,3) = 5
              ENDIF
              obj(nobj)%npts(3) = 4
              obj(nobj)%ipts(1,3) = 3
              obj(nobj)%ipts(2,3) = 2
              obj(nobj)%ipts(3,3) = 6
              obj(nobj)%ipts(4,3) = 7
c...          Surface 4:
              IF (iy.EQ.nybin) THEN   
                obj(nobj)%tsur(4) = SP_GRID_BOUNDARY
                obj(nobj)%nmap(4) = 1
                obj(nobj)%imap(1,4) = nobj
                obj(nobj)%isur(1,4) = 4
              ELSE
                obj(nobj)%tsur(4) = SP_GRID_SURFACE
                obj(nobj)%nmap(4) = 1
                obj(nobj)%imap(1,4) = nobj - nxbin
                obj(nobj)%isur(1,4) = 2
              ENDIF
              obj(nobj)%npts(4) = 4
              obj(nobj)%ipts(1,4) = 4    ! If this is right, then the rotation is backwards! 
              obj(nobj)%ipts(2,4) = 3
              obj(nobj)%ipts(3,4) = 7
              obj(nobj)%ipts(4,4) = 8
c...          Surface 5:
              IF (ix.EQ.nxbin) THEN
                obj(nobj)%tsur(5) = SP_GRID_BOUNDARY
                obj(nobj)%nmap(5) = 1
                obj(nobj)%imap(1,5) = nobj
                obj(nobj)%isur(1,5) = 5
              ELSE
                obj(nobj)%tsur(5) = SP_GRID_SURFACE
                obj(nobj)%nmap(5) = 1
                obj(nobj)%imap(1,5) = nobj + 1
                obj(nobj)%isur(1,5) = 3
              ENDIF
              obj(nobj)%npts(5) = 4
              obj(nobj)%ipts(1,5) = 1
              obj(nobj)%ipts(2,5) = 4
              obj(nobj)%ipts(3,5) = 8
              obj(nobj)%ipts(4,5) = 5
c...          Surface 6:
              IF (isector.EQ.iend) THEN
                IF (iend-istart+1.EQ.nsector) THEN   
c...              Full torus:
                  obj(nobj)%tsur(6) = SP_GRID_SURFACE
                  obj(nobj)%nmap(6) = 1
                  obj(nobj)%imap(1,6) = nobj - ncells * (nsector - 1)
                  obj(nobj)%isur(1,6) = 1
                ELSE
                  obj(nobj)%tsur(6) = SP_GRID_BOUNDARY
                  obj(nobj)%nmap(6) = 1
                  obj(nobj)%imap(1,6) = nobj
                  obj(nobj)%isur(1,6) = 6
                ENDIF
              ELSE
                obj(nobj)%tsur(6) = SP_GRID_SURFACE
                obj(nobj)%nmap(6) = 1
                obj(nobj)%imap(1,6) = nobj + ncells
                obj(nobj)%isur(1,6) = 1
              ENDIF
              obj(nobj)%npts(6) = 4
              obj(nobj)%ipts(1,6) = 8
              obj(nobj)%ipts(2,6) = 7
              obj(nobj)%ipts(3,6) = 6
              obj(nobj)%ipts(4,6) = 5
            ENDDO
          ENDDO
 
        ENDDO

      ELSE
        CALL ER('BuildInversionMesh','Unknown option',*99)
      ENDIF

c...  Local mesh refinement adjustment based on mask, from a previous 
c     iteration:


c...  Generalized cell mapping method:

c      WRITE(0,*) '???',obj(8)%gsur(1)


 98   RETURN
 99   STOP
      END

