c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE LoadEmissionDatabase(dummy,ni,nr)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

 
      CHARACTER dummy*(*)
      INTEGER   ni,nr

      INTEGER   i1,idum1
      REAL      rdum1
      CHARACTER cdum1*256

      SELECTCASE (opt%int_database(opt%int_num))
        CASE (1,3,4)
          READ(dummy,*) cdum1,(idum1,i1=1,ni),
     .                        (rdum1,i1=1,nr),idum1,
     .                  opt%int_line(opt%int_num)
        CASE (2)
          READ(dummy,*) cdum1,(idum1,i1=1,ni),
     .                        (rdum1,i1=1,nr),idum1,
     .                  opt%int_adasid(opt%int_num),
     .                  opt%int_adasyr(opt%int_num),
     .                  opt%int_adasex(opt%int_num),
     .                  opt%int_isele (opt%int_num),
     .                  opt%int_iselr (opt%int_num),
     .                  opt%int_iselx (opt%int_num),
     .                  opt%int_iseld (opt%int_num)
        CASE DEFAULT
          STOP 'SORRY, NO USER OPTIONS YET FOR DATABASE'
      ENDSELECT

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      LOGICAL FUNCTION GetLine(fp,buffer,mode)
      IMPLICIT none

      INTEGER  , INTENT(IN)  :: fp,mode
      CHARACTER, INTENT(OUT) :: buffer*(*)

      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

      INTEGER i,n

      GetLine = .TRUE. 

      DO WHILE(.TRUE.)
        READ(fp,'(A)',END=10) buffer

c...    Remove leading spaces:
        n = LEN_TRIM(buffer)
        DO i = 1, n
          IF (buffer(i:i).NE.' ') EXIT
        ENDDO 
        IF (i.LT.n) buffer(1:n) = buffer(i:n+i)      

c...    Remove portion of line after comment charcter:
        n = LEN_TRIM(buffer)
        DO i = 1, n
          IF (buffer(i:i).EQ.'$'.OR.buffer(i:i).EQ.'*') EXIT
        ENDDO 
        buffer(i:n) = ' '

        IF (LEN_TRIM(buffer).EQ.0) THEN
c...      Comment or blank line, so continue:
        ELSEIF (buffer(1:1).EQ.'{') THEN
          IF (mode.EQ.WITH_TAG) THEN
c...        Tag line: 
            EXIT
          ELSE
c...        Tag line was not requested so backup the file position:
            GetLine = .FALSE.
            BACKSPACE(fp)
            EXIT
          ENDIF
        ELSE
          IF (mode.EQ.NO_TAG) THEN
c...        Data line found, as requested:
            EXIT
          ELSE
c...        Data line found but looking for next tag, so keep going:
          ENDIF
        ENDIF
      ENDDO
c      WRITE(0,*) 'BUFFER:',buffer(1:50),GetLine
      RETURN

 10   GetLine = .FALSE.
      RETURN
      END
c
c ======================================================================
c
      SUBROUTINE LoadOptions985_New(opt,mode,status)
      USE mod_out985
      USE mod_out989
      IMPLICIT none

      INTEGER, INTENT(IN)  :: mode
      INTEGER, INTENT(OUT) :: status
      TYPE(type_options985), INTENT(OUT) :: opt     

      LOGICAL GetLine

      INTEGER       i1,i2,fp,idum1,idum2,i,n,inext
      LOGICAL       load_detector
      CHARACTER     cdum1*128,dummy*1024,buffer*1024
      CHARACTER*256 buffer_array(100)

      TYPE(type_header) :: header
      REAL*8 :: image(1000,1000)
      INTEGER ix,iy
      TYPE(type_options989) :: opt989      

      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

      WRITE(0,*) '  MAIN985 New...'

      status = 0

      fp = 5

      load_detector = .FALSE.
      IF (mode.EQ.DETECTOR_ONLY) load_detector = .TRUE.

c...  Scan input file:
      DO WHILE(GetLine(fp,buffer,WITH_TAG))
c        READ(fp,'(A256)',END=10) buffer      

        n = LEN_TRIM(buffer)

c        WRITE(0,*) 'BUFFER:',TRIM(buffer)

c        IF (buffer(2:2).NE.'{') CYCLE

c...    Remove the portion of the data line that comes after a comment
c       character:
c        DO j = 1, LEN_TRIM(buffer)                               ! THIS SHOULD ALL BE REPLACED 
c          IF (buffer(j:j).EQ.'$'.OR.buffer(j:j).EQ.'*') EXIT     ! BY A GETLINE CALL IN THE LOOP I THINK...
c        ENDDO 
c        buffer(j:LEN_TRIM(buffer)) = ' '

c...    Isolate tag string:
        DO i = 2, LEN_TRIM(buffer)
          IF (buffer(i:i).EQ.'}') EXIT
        ENDDO

        SELECTCASE (buffer(2:i-1))
c         --------------------------------------------------------------
          CASE('RAY TRACE') 
c         --------------------------------------------------------------
          CASE('GEOMETRY')
            IF (mode.NE.ALL_OPTIONS) CYCLE
            READ(buffer(i+1:n),*) idum1
c            WRITE(0,*) 'BUFFER:-',buffer(i+1:n),idum1
            WRITE(0,*) 'LOADING GEOMETRY:',idum1
            IF (idum1.NE.0) THEN
              DO WHILE(GetLine(fp,buffer,NO_TAG))
                READ(buffer,*) cdum1,idum1
c                WRITE(0,*) 'BUFFER:-',TRIM(buffer)
                IF (idum1.EQ.0) CYCLE
c                WRITE(0,*) 'LOADING GEOMETRY:',idum1
                opt%obj_num = opt%obj_num + 1
                SELECTCASE (idum1)
                  CASE (1)  ! Data from magnetic/standard OSM fluid grid
                    READ(buffer,*) cdum1,
     .                opt%obj_type  (opt%obj_num),
     .                opt%obj_option(opt%obj_num),
     .                opt%obj_colour(opt%obj_num),
     .                opt%obj_reflec(opt%obj_num),
     .                opt%obj_fudge (opt%obj_num),
     .                opt%obj_factor(opt%obj_num)
                  CASE (2)  ! Data from Eirene triangle grid
                    READ(buffer,*) cdum1,
     .                opt%obj_type  (opt%obj_num),
     .                opt%obj_option(opt%obj_num),
     .                opt%obj_colour(opt%obj_num),
     .                opt%obj_reflec(opt%obj_num),
     .                opt%obj_fudge (opt%obj_num),
     .                opt%obj_factor(opt%obj_num),
     .                opt%obj_n     (opt%obj_num,1),
     .                opt%obj_n     (opt%obj_num,2)
                  CASE (3)  ! Inversion mesh
                    READ(buffer,*) cdum1,
     .                opt%obj_type  (opt%obj_num),
     .                opt%obj_option(opt%obj_num),
     .                opt%obj_colour(opt%obj_num),
     .                opt%obj_reflec(opt%obj_num),
     .                opt%obj_fudge (opt%obj_num),
     .                opt%obj_factor(opt%obj_num),
     .                opt%obj_n     (opt%obj_num,1),
     .                opt%obj_n     (opt%obj_num,2),
     .                opt%obj_r     (opt%obj_num,1),
     .                opt%obj_z     (opt%obj_num,1),
     .                opt%obj_r     (opt%obj_num,2),
     .                opt%obj_z     (opt%obj_num,2)
                  CASE (4)  ! Raw triangle data from drawings office
                    READ(buffer,*) cdum1,
     .                opt%obj_type  (opt%obj_num),
     .                opt%obj_option(opt%obj_num),
     .                opt%obj_colour(opt%obj_num),
     .                opt%obj_reflec(opt%obj_num),
     .                opt%obj_fudge (opt%obj_num),
     .                opt%obj_factor(opt%obj_num),
     .                opt%obj_n     (opt%obj_num,1),
     .                opt%obj_n     (opt%obj_num,2),
     .                opt%obj_orientation(opt%obj_num),
     .                opt%obj_scale      (opt%obj_num),
     .                opt%obj_yangle     (opt%obj_num),
     .                opt%obj_fname      (opt%obj_num)
                  CASE (5)  ! File with a list of line segments:
                    READ(buffer,*) cdum1,
     .                opt%obj_type  (opt%obj_num),
     .                opt%obj_option(opt%obj_num),
     .                opt%obj_colour(opt%obj_num),
     .                opt%obj_reflec(opt%obj_num),
     .                opt%obj_fudge (opt%obj_num),
     .                opt%obj_factor(opt%obj_num),
     .                opt%obj_n     (opt%obj_num,1),
     .                opt%obj_n     (opt%obj_num,2),
     .                opt%obj_fname (opt%obj_num)
                  CASE (6)  ! Data from tetrahedron grid:
                    READ(buffer,*) cdum1,
     .                opt%obj_type  (opt%obj_num),
     .                opt%obj_option(opt%obj_num),
     .                opt%obj_colour(opt%obj_num),
     .                opt%obj_reflec(opt%obj_num),
     .                opt%obj_fudge (opt%obj_num),
     .                opt%obj_factor(opt%obj_num),
     .                opt%obj_n     (opt%obj_num,1),
     .                opt%obj_n     (opt%obj_num,2),
     .                opt%obj_fname (opt%obj_num)
                  CASE (7)  ! Reconstructed image
                    READ(buffer,*) cdum1,
     .                opt%obj_type  (opt%obj_num),
     .                opt%obj_option(opt%obj_num),
     .                opt%obj_colour(opt%obj_num),
     .                opt%obj_reflec(opt%obj_num),
     .                opt%obj_fudge (opt%obj_num),
     .                opt%obj_factor(opt%obj_num),
     .                opt%obj_n     (opt%obj_num,1),
     .                opt%obj_n     (opt%obj_num,2),
     .                opt%obj_fname (opt%obj_num)
                  CASE (8)  ! ITER first wall panel
                    READ(buffer,*) cdum1,
     .                opt%obj_type  (opt%obj_num),
     .                opt%obj_option(opt%obj_num),
     .                opt%obj_colour(opt%obj_num),
     .                opt%obj_reflec(opt%obj_num),
     .                opt%obj_fudge (opt%obj_num),
     .                opt%obj_factor(opt%obj_num),
     .                opt%obj_n     (opt%obj_num,1),
     .                opt%obj_n     (opt%obj_num,2),
     .                opt%obj_sym   (opt%obj_num),
     .                opt%obj_fname (opt%obj_num)
                  CASE DEFAULT
                    CALL User_LoadVesselGeometry(opt,idum1,buffer)
                ENDSELECT
              ENDDO
            ENDIF
c         --------------------------------------------------------------
          CASE('REFLECTIONS')
            IF (mode.NE.ALL_OPTIONS) CYCLE
            READ(buffer(i+1:n),*) idum1
            WRITE(0,*) 'LOADING REFLECTIONS:',idum1
            IF (idum1.NE.0) THEN
              DO WHILE(GetLine(fp,buffer,NO_TAG))
                READ(buffer,*) cdum1,idum1
                IF (idum1.EQ.0) CYCLE
                opt%ref_num = opt%ref_num + 1
                SELECTCASE (idum1)
                  CASE (1)  ! Specular 
                    READ(buffer,*) cdum1,
     .                opt%ref_model (opt%ref_num),  ! 1=specular, 2=diffuse
     .                opt%ref_wlgth (opt%ref_num),  ! Variation of reflectivity
     .                opt%ref_k     (opt%ref_num),  ! Drop in signal for uniform reflectivity option
     .                opt%ref_ow    (opt%ref_num),  ! Not in use...
     .                opt%ref_pw    (opt%ref_num),  ! Not in use...
     .                opt%ref_n     (opt%ref_num),  ! Some kind of exponent
     .                opt%ref_cutoff(opt%ref_num),  ! Don't include contibutions below a certain relative intensity
     .                opt%ref_otheta(opt%ref_num),  ! THETA angle distribution option (0=none, 1=linear)
     .                opt%ref_dtheta(opt%ref_num),  ! THETA angle step size
     .                opt%ref_ophi  (opt%ref_num),  ! PHI angle distribution option (0=none, 1=linear)
     .                opt%ref_dphi  (opt%ref_num)   ! PHI angle step size
                  CASE (2)  ! Diffuse
                    READ(buffer,*) cdum1,
     .                opt%ref_model (opt%ref_num),
     .                opt%ref_wlgth (opt%ref_num),
     .                opt%ref_k     (opt%ref_num),
     .                opt%ref_otheta(opt%ref_num),
     .                opt%ref_dtheta(opt%ref_num),
     .                opt%ref_ophi  (opt%ref_num),
     .                opt%ref_dphi  (opt%ref_num)
                  CASE (3)
                  CASE DEFAULT
                ENDSELECT
              ENDDO
            ENDIF
c         --------------------------------------------------------------
          CASE('INTEGRATION')
            IF (mode.NE.ALL_OPTIONS) CYCLE
            READ(buffer(i+1:n),*) idum1
            IF (idum1.NE.0) THEN
              DO WHILE(GetLine(fp,buffer,NO_TAG))
                READ(buffer,*) cdum1,idum1
                IF (idum1.EQ.0) CYCLE
                WRITE(0,*) 'LOADING INTEGRATION:',idum1
                opt%int_num = opt%int_num + 1
                SELECTCASE (idum1)
                  CASE (1)  ! Straight-up line integral:
                    READ(buffer,*) cdum1,opt%int_type    (opt%int_num),                
     .                                   opt%int_colour  (opt%int_num),
     .                                   opt%int_z       (opt%int_num),
     .                                   opt%int_a       (opt%int_num),
     .                                   opt%int_charge  (opt%int_num),
     .                                   opt%int_index   (opt%int_num),
     .                                   opt%int_database(opt%int_num)
c...                Data source:
                    CALL LoadEmissionDatabase(buffer,6,0)
                  CASE (2)  ! Line shape integral:
                    READ(buffer,*) cdum1,opt%int_type    (opt%int_num),                
     .                                   opt%int_colour  (opt%int_num),
     .                                   opt%int_z       (opt%int_num),
     .                                   opt%int_a       (opt%int_num),
     .                                   opt%int_charge  (opt%int_num),
     .                                   opt%int_index   (opt%int_num),
     .                                   opt%int_shape   (opt%int_num),
     .                                   opt%int_instr   (opt%int_num),
     .                                   opt%int_width   (opt%int_num),
     .                                   opt%int_database(opt%int_num)
c...                Data source:
                    CALL LoadEmissionDatabase(buffer,7,2)
                  CASE (3)  ! Line-of-sight weighted average:
                    READ(buffer,*) cdum1,opt%int_type    (opt%int_num),                
     .                                   opt%int_colour  (opt%int_num),
     .                                   opt%int_z       (opt%int_num),
     .                                   opt%int_a       (opt%int_num),
     .                                   opt%int_charge  (opt%int_num),
     .                                   opt%int_index   (opt%int_num),
     .                                   opt%int_average (opt%int_num),
     .                                   opt%int_database(opt%int_num)
c...                Data source:
                    CALL LoadEmissionDatabase(buffer,7,0)
                  CASE DEFAULT
                    CALL ER('LoadOptions985_New','Unrecognized '//
     .                      'integration option',*99)
                ENDSELECT
              ENDDO
            ENDIF
c         --------------------------------------------------------------
          CASE('DETECTOR')
            IF (mode.NE.DETECTOR_ONLY.OR..NOT.load_detector) THEN
              BACKSPACE(fp)
              EXIT
            ENDIF
            READ(buffer(i+1:n),*) idum1
            WRITE(0,*) 'LOADING DETECTOR:',idum1
            SELECTCASE(idum1)
              CASE(0)
              CASE(1:2)
                load_detector = .FALSE.
                opt%ccd = idum1
                READ(fp,*) cdum1,opt%focallength
                READ(fp,*) cdum1,opt%distortion
                READ(fp,*) cdum1,opt%cen(1:3)
                READ(fp,*) cdum1,opt%roll,opt%tilt,opt%swing
                READ(fp,*) cdum1,opt%width(1:2)
                READ(fp,*) cdum1,opt%angle(1:2)  ! Maybe only need to load the first entry
                READ(fp,*) cdum1,opt%nxbin,opt%nybin
                READ(fp,*) cdum1,opt%sa_nxbin,opt%sa_nybin,opt%sa_opt,
     .                           opt%sa_par1 ,opt%sa_par2
                READ(fp,*) cdum1,opt%fmap
              CASE(3:4)
                load_detector = .FALSE.
                opt%ccd = idum1 - 2
                 
              CASE DEFAULT
                CALL ER('LoadOptions985_New','Unknown DETECTOR '//
     .                  'option',*99)          
            ENDSELECT
c         --------------------------------------------------------------
          CASE('DETECTOR MASK')
            IF (mode.NE.DETECTOR_ONLY.OR.load_detector) CYCLE            
            READ(buffer(i+1:n),*) idum1
            WRITE(0,*) 'LOADING MASK:',idum1
            IF (idum1.NE.0) THEN
              DO WHILE(GetLine(fp,buffer,NO_TAG))
                READ(buffer,*) cdum1,idum1
                IF (idum1.EQ.0) CYCLE
                opt%mask_num = opt%mask_num + 1
                SELECTCASE (idum1)
                  CASE (-1,1)
c...                Load image and look for pure/artificial black (pixel value=0):
                    READ(buffer,*) cdum1,opt%mask_opt(opt%mask_num),
     .                                   opt989%imagetype,
     .                                   opt989%fimage
                    WRITE(0,*) 'LOADING IMAGE:',TRIM(opt989%fimage)
                    CALL LoadImage(opt989,image,header,1000,1000)  ! find out image size... as elsewhere... 
                    opt%img_opt = 1
                    opt%img_nxbin = opt989%nxbin
                    opt%img_nybin = opt989%nybin
                    DO ix = 1, opt989%nxbin
                      DO iy = 1, opt989%nybin
                        opt%img_image(ix,iy) = image(ix,iy)  ! Was using a pointer but SUN compiler complained...
                      ENDDO
                    ENDDO
                    WRITE(0,*) 'nxbin,nybin:',
     .                         opt%img_nxbin,opt%img_nybin
                    WRITE(0,*) '...',opt%mask_num,
     .                         opt%mask_opt(opt%mask_num)
     .                         
                  CASE (2)  ! Invert mask
                    READ(buffer,*) cdum1,opt%mask_opt(opt%mask_num)
                  CASE (3)
                  CASE (4)  ! Line
                    READ(buffer,*) cdum1,opt%mask_opt(opt%mask_num),
     .                                   opt%mask_vtx(opt%mask_num,1:4)
                  CASE (5)  ! Box
                    READ(buffer,*) cdum1,opt%mask_opt(opt%mask_num),
     .                                   opt%mask_vtx(opt%mask_num,1:4)
                  CASE DEFAULT
                    CALL ER('LoadOptions985_New','Unrecognized '//
     .                      'detector mask type',*99)
                ENDSELECT
              ENDDO
            ENDIF
c         --------------------------------------------------------------
          CASE('RIBBON')
            READ(buffer(i+1:n),*) idum1
            WRITE(0,*) 'LOADING RIBBON:',idum1
            IF (idum1.LE.0) CYCLE
            DO WHILE(GetLine(fp,buffer,NO_TAG))
              READ(buffer,*) cdum1,idum1
              IF (idum1.EQ.0) CYCLE
              CALL SplitBuffer(buffer,buffer_array) 
              opt%rib_n = opt%rib_n + 1
              i1 = opt%rib_n
              SELECTCASE (idum1)
                CASE (1)  ! 
                  opt%rib_option(i1) = idum1
                  READ(buffer_array(3),*) opt%rib_nrad (  i1)
                  READ(buffer_array(4),*) opt%rib_nphi (  i1)
                  READ(buffer_array(5),*) opt%rib_phi  (1,i1)
                  READ(buffer_array(6),*) opt%rib_phi  (2,i1)
                  READ(buffer_array(7),*) opt%rib_scale(  i1)
                  READ(buffer_array(8),*) opt%rib_trace(  i1)
                  SELECTCASE (opt%rib_trace(i1))
                    CASE (1) 
                      READ(buffer_array(9 ),*) opt%rib_dphi(i1)
                      READ(buffer_array(10),*) opt%rib_r   (1,i1)
                      READ(buffer_array(11),*) opt%rib_z   (1,i1)
                      READ(buffer_array(12),*) opt%rib_r   (2,i1)
                      READ(buffer_array(13),*) opt%rib_z   (2,i1)
                      inext = 14
                    CASE (2) 
                      opt%rib_tfile(i1) = TRIM(buffer_array(9))
                      inext = 10
                    CASE DEFAULT
                      CALL ER('LoadOptions985_New','Unrecognized '//
     .                        'trace option',*99)
                  ENDSELECT
                  READ(buffer_array(inext),*) opt%rib_limit(i1)
                  SELECTCASE (opt%rib_limit(i1))
                    CASE (0) 
                      opt%rib_tag(i1) = TRIM(buffer_array(inext+1))
                    CASE DEFAULT
                      CALL ER('LoadOptions985_New','Unrecognized '//
     .                        'ribbon limit',*99)
                  ENDSELECT
                CASE DEFAULT
                  CALL ER('LoadOptions985_New','Unrecognized '//
     .                    'ribbon grid option',*99)
              ENDSELECT
            ENDDO
c         --------------------------------------------------------------
          CASE('PLOTS')
            IF (mode.NE.ALL_OPTIONS) CYCLE
            IF (opt%nplots.EQ.-1) opt%nplots = 0
            READ(buffer(i+1:n),*) idum1
            IF (idum1.NE.0) THEN
              DO WHILE(GetLine(fp,buffer,NO_TAG))
                opt%nplots = opt%nplots + 1
                opt%plots(opt%nplots) = buffer(1:LEN_TRIM(buffer))  ! Still don't know why TRIM() doesn't work...
c                WRITE(0,*) opt%nplots,buffer(1:LEN_TRIM(buffer))
              ENDDO
            ENDIF
c         --------------------------------------------------------------
          CASE('END')
            status = -1
            EXIT
c         --------------------------------------------------------------
          CASE DEFAULT
            CALL ER('LoadOptions985_New','Unrecognized tag',*99)
        ENDSELECT

      ENDDO


      RETURN
 99   WRITE(0,*) '  BUFFER: >'//buffer(1:50)//'<'
      STOP
      END
c
c ======================================================================
c
c subroutine: LoadOptions985
c
      SUBROUTINE LoadOptions985(opt)
      USE mod_out985
      IMPLICIT none

      TYPE(type_options985) :: opt     

      INTEGER i1,i2,fp,idum1,idum2
      CHARACTER cdum1*128,dummy*1024

      WRITE(0,*) '  MAIN985'

      READ(5,'(A256)') dummy
      WRITE(0,*) '>>'//dummy(8:15)//'<<'
      IF (dummy(8:15).EQ.'3D MODEL'.OR.dummy(8:15).EQ.'3D Model'.OR.
     .    dummy(8:15).EQ.'3d model'.OR.dummy(8:15).EQ.'3D model') THEN
        WRITE(0,*) 'LOADING MODEL DATA FROM FILE'
        BACKSPACE 5
c        READ(5,*) cdum1,opt%ob_model
c        READ(5,*) cdum1,
c     .            opt%ob_nsector,opt%ob_angle_start,opt%ob_angle_end,
c     .            opt%ob_yrotation
c        READ(5,*) cdum1,
c     .            opt%ob_stdgrd,opt%ob_stdgrd_colour,
c     .            opt%ob_stdgrd_reflec
c        READ(5,*) cdum1,
c     .            opt%ob_trigrd,opt%ob_trigrd_colour,
c     .            opt%ob_trigrd_reflec
c        READ(5,*) cdum1,
c     .    opt%ob_invgrd,opt%ob_invgrd_colour,opt%ob_invgrd_idum1,
c     .    opt%ob_invgrd_xcen  ,opt%ob_invgrd_ycen,
c     .    opt%ob_invgrd_xwidth,opt%ob_invgrd_ywidth,
c     .    opt%ob_invgrd_nxbin ,opt%ob_invgrd_nybin
c        READ(5,*) cdum1,opt%ob_wall,opt%ob_wall_colour,
c     .            opt%ob_wall_reflec
c        READ(5,*) cdum1,opt%ob_targ,opt%ob_targ_colour,
c     .            opt%ob_targ_reflec

c        opt%ob_raw_num = 0
c        DO i1 = 1, 10
c          READ(5,'(A256)') dummy
cc          WRITE(0,*)   '>???>'//dummy(8:19)//'<<'
c          IF (dummy(8:13).EQ.'  raw '.OR.dummy(8:13).EQ.'  RAW '.OR.
c     .        dummy(8:13).EQ.'  Raw '.OR.dummy(8:13).EQ.'      ') THEN
c            READ(dummy,*) cdum1,idum1
c            IF (idum1.NE.0) THEN
c              opt%ob_raw_num = opt%ob_raw_num + 1
c              IF (opt%ob_raw_num.GT.10) 
c     .          CALL ER('Main985','Whao! to much raw',*99)
c              READ(dummy,*) cdum1,idum1,
c     .          opt%ob_raw_colour     (opt%ob_raw_num),
c     .          opt%ob_raw_reflec     (opt%ob_raw_num),
c     .          opt%ob_raw_ind        (opt%ob_raw_num,1),
c     .          opt%ob_raw_ind        (opt%ob_raw_num,2),
c     .          opt%ob_raw_material   (opt%ob_raw_num),
c     .          opt%ob_raw_orientation(opt%ob_raw_num),  ! Also need reflection about midplane, and rotation... 
c     .          opt%ob_raw_scale      (opt%ob_raw_num), 
c     .          opt%ob_raw_fname      (opt%ob_raw_num)
cc...          Check for custom reflection model:
c            ENDIF
c          ELSE
c            BACKSPACE 5
c            EXIT
c          ENDIF
c        ENDDO

c        READ(5,*) cdum1,opt%ob_tube,opt%ob_tube_colour,opt%ob_tube_idum1
c        READ(5,*) cdum1,opt%ob_line,opt%ob_line_colour,opt%ob_line_idum1
c        DO i1 = 1, 10
c          READ(5,'(A256)') dummy
c          IF (dummy(8:13).EQ.'  user'.OR.dummy(8:13).EQ.'  user'.OR.
c     .        dummy(8:13).EQ.'  user') THEN          
c            READ(dummy,*) cdum1,opt%ob_user(i1),
c     .                    opt%ob_user_colour(i1),
c     .                    opt%ob_user_reflec(i1)
c            WRITE(0,*) 'USER FOUND:',i1
c          ELSE
c            BACKSPACE 5
c            EXIT
c          ENDIF
c        ENDDO
c...    Check input to make sure requested models are conistent:

      ELSE
        BACKSPACE 5
c...    Defaults:
c        opt%ob_model = 1
c        opt%ob_nsector = 48    
c        opt%ob_angle_start = 0.0
c        opt%ob_angle_end = 360.0
c        opt%ob_yrotation = 0.0
c        opt%ob_stdgrd = 1
c        opt%ob_stdgrd_colour = 3
c        opt%ob_trigrd = 0
c        opt%ob_invgrd = 0
c        opt%ob_wall = 1
c        opt%ob_wall_colour = 1
c        opt%ob_targ = 0
c        opt%ob_tube = 0
c        opt%ob_line = 0
c        opt%ob_user = 0
      ENDIF


      opt%obj_num = 0
      READ(5,'(A256)') dummy
      WRITE(0,*) '>>'//dummy(8:15)//'<<'
      IF (dummy(8:14).EQ.'GENERIC'.OR.dummy(8:15).EQ.'Generic'.OR.
     .    dummy(8:14).EQ.'generic') THEN
        WRITE(0,*) 'LOADING GENERIC 3D MODEL DATA'
        READ(dummy,*) cdum1,idum1
        DO i1 = 1, 20
          READ(5,'(A256)') dummy
          IF (dummy(8:9).EQ.'  ') THEN
            IF (idum1.NE.0) THEN
              READ(dummy,*) cdum1,idum2
              IF (idum2.NE.0) opt%obj_num = opt%obj_num + 1
              SELECTCASE (idum2)
                CASE (0)
                CASE (1)  ! Data from magnetic/standard OSM grid
                  READ(dummy,*) cdum1,
     .              opt%obj_type  (opt%obj_num),
     .              opt%obj_option(opt%obj_num),
     .              opt%obj_colour(opt%obj_num),
     .              opt%obj_reflec(opt%obj_num),
     .              opt%obj_fudge (opt%obj_num),
     .              opt%obj_factor(opt%obj_num)
                CASE (2)  ! Data from EIRENE triangle grid
                  READ(dummy,*) cdum1,
     .              opt%obj_type  (opt%obj_num),
     .              opt%obj_option(opt%obj_num),
     .              opt%obj_colour(opt%obj_num),
     .              opt%obj_reflec(opt%obj_num),
     .              opt%obj_fudge (opt%obj_num),
     .              opt%obj_factor(opt%obj_num),
     .              opt%obj_n     (opt%obj_num,1),
     .              opt%obj_n     (opt%obj_num,2)
                CASE (3)  ! Inversion mesh
                  READ(dummy,*) cdum1,
     .              opt%obj_type  (opt%obj_num),
     .              opt%obj_option(opt%obj_num),
     .              opt%obj_colour(opt%obj_num),
     .              opt%obj_reflec(opt%obj_num),
     .              opt%obj_fudge (opt%obj_num),
     .              opt%obj_factor(opt%obj_num),
     .              opt%obj_n     (opt%obj_num,1),
     .              opt%obj_n     (opt%obj_num,2),
     .              opt%obj_r     (opt%obj_num,1),
     .              opt%obj_z     (opt%obj_num,1),
     .              opt%obj_r     (opt%obj_num,2),
     .              opt%obj_z     (opt%obj_num,2)
                CASE (4)  ! Raw triangle data from drawings office
                  READ(dummy,*) cdum1,
     .              opt%obj_type  (opt%obj_num),
     .              opt%obj_option(opt%obj_num),
     .              opt%obj_colour(opt%obj_num),
     .              opt%obj_reflec(opt%obj_num),
     .              opt%obj_fudge (opt%obj_num),
     .              opt%obj_factor(opt%obj_num),
     .              opt%obj_n     (opt%obj_num,1),
     .              opt%obj_n     (opt%obj_num,2),
     .              opt%obj_orientation(opt%obj_num),
     .              opt%obj_scale      (opt%obj_num),
     .              opt%obj_yangle     (opt%obj_num),
     .              opt%obj_fname      (opt%obj_num)
                CASE (5)  ! File with a list of line segments:
                  READ(dummy,*) cdum1,
     .              opt%obj_type  (opt%obj_num),
     .              opt%obj_option(opt%obj_num),
     .              opt%obj_colour(opt%obj_num),
     .              opt%obj_reflec(opt%obj_num),
     .              opt%obj_fudge (opt%obj_num),
     .              opt%obj_factor(opt%obj_num),
     .              opt%obj_n     (opt%obj_num,1),
     .              opt%obj_n     (opt%obj_num,2),
     .              opt%obj_fname (opt%obj_num)
                CASE (6)  ! Data from tetrahedron grid:
                  READ(dummy,*) cdum1,
     .              opt%obj_type  (opt%obj_num),
     .              opt%obj_option(opt%obj_num),
     .              opt%obj_colour(opt%obj_num),
     .              opt%obj_reflec(opt%obj_num),
     .              opt%obj_fudge (opt%obj_num),
     .              opt%obj_factor(opt%obj_num)
                CASE DEFAULT
                  CALL User_LoadVesselGeometry(opt,idum2,dummy)
              ENDSELECT
            ENDIF
          ELSE
            BACKSPACE 5
            EXIT
          ENDIF
        ENDDO
        WRITE(0,*) '  OBJ_NUM=',opt%obj_num
      ELSE
        BACKSPACE 5
c...    Defaults:
      ENDIF


      opt%ref_num = 0
      READ(5,'(A256)') dummy
      WRITE(0,*) '>>'//dummy(8:15)//'<<'
      IF (dummy(8:15).EQ.'REFLECTI'.OR.dummy(8:15).EQ.'Reflecti'.OR.
     .    dummy(8:15).EQ.'reflecti') THEN
        WRITE(0,*) 'LOADING REFLECTION MODEL DATA'
        READ(dummy,*) cdum1,idum1
        DO i1 = 1, 10
          READ(5,'(A256)') dummy
c          WRITE(0,*)   '>REF>'//dummy(8:19)//'<<'
          IF (dummy(8:9).EQ.'  ') THEN
            IF (idum1.NE.0) THEN
              READ(dummy,*) cdum1,idum2
              IF (idum2.NE.0) opt%ref_num = opt%ref_num + 1
              SELECTCASE (idum2)
                CASE (1)  ! Specular reflection model
                  READ(dummy,*) cdum1,
     .              opt%ref_model (opt%ref_num),
     .              opt%ref_wlgth (opt%ref_num),
     .              opt%ref_k     (opt%ref_num),
     .              opt%ref_ow    (opt%ref_num),
     .              opt%ref_pw    (opt%ref_num),
     .              opt%ref_n     (opt%ref_num),
     .              opt%ref_cutoff(opt%ref_num),
     .              opt%ref_otheta(opt%ref_num),
     .              opt%ref_dtheta(opt%ref_num),
     .              opt%ref_ophi  (opt%ref_num),
     .              opt%ref_dphi  (opt%ref_num)
                CASE (2)  ! Diffuse
                  READ(dummy,*) cdum1,
     .              opt%ref_model (opt%ref_num),
     .              opt%ref_wlgth (opt%ref_num),
     .              opt%ref_k     (opt%ref_num),
     .              opt%ref_otheta(opt%ref_num),
     .              opt%ref_dtheta(opt%ref_num),
     .              opt%ref_ophi  (opt%ref_num),
     .              opt%ref_dphi  (opt%ref_num)
                CASE (3)
                CASE DEFAULT
              ENDSELECT
            ENDIF
          ELSE
            BACKSPACE 5
            EXIT
          ENDIF
        ENDDO
        WRITE(0,*) '  REF_NUM=',opt%ref_num
        DO i1 = 1, opt%ref_num
          WRITE(0,'(2I4,F7.2,4(I4,F7.2,2X))')
     .      opt%ref_model (opt%ref_num),
     .      opt%ref_wlgth (opt%ref_num),
     .      opt%ref_k     (opt%ref_num),
     .      opt%ref_ow    (opt%ref_num),
     .      opt%ref_pw    (opt%ref_num),
     .      opt%ref_n     (opt%ref_num),
     .      opt%ref_cutoff(opt%ref_num),
     .      opt%ref_otheta(opt%ref_num),
     .      opt%ref_dtheta(opt%ref_num),
     .      opt%ref_ophi  (opt%ref_num),
     .      opt%ref_dphi  (opt%ref_num)
        ENDDO
      ELSE
        BACKSPACE 5
c...    Defaults:
      ENDIF



      opt%int_num = 0
      READ(5,'(A256)') dummy
      WRITE(0,*) '>>'//dummy(8:15)//'<<'
      IF (dummy(8:15).EQ.'INTEGRAT'.OR.dummy(8:15).EQ.'Integrat'.OR.
     .    dummy(8:15).EQ.'integrat') THEN
        WRITE(0,*) 'LOADING INTEGRATION DATA'
        DO i1 = 1, 10
          READ(5,'(A256)') dummy
c          WRITE(0,*)   '>REF>'//dummy(8:19)//'<<'
          IF (dummy(8:9).EQ.'  ') THEN
            READ(dummy,*) cdum1,idum1
            opt%int_num = opt%int_num + 1
            SELECTCASE (idum1)
              CASE (0)
c...            Turned off:
                opt%int_num = opt%int_num - 1
              CASE (1)  
c...            Straight-up line integral:
                READ(dummy,*) cdum1,opt%int_type    (opt%int_num),                
     .                              opt%int_colour  (opt%int_num),
     .                              opt%int_z       (opt%int_num),
     .                              opt%int_a       (opt%int_num),
     .                              opt%int_charge  (opt%int_num),
     .                              opt%int_index   (opt%int_num),
     .                              opt%int_database(opt%int_num)
c...            Data source:
                CALL LoadEmissionDatabase(dummy,6,0)

              CASE (2)
c...            Line shape integral:
c...            Straight-up line integral:
                READ(dummy,*) cdum1,opt%int_type    (opt%int_num),                
     .                              opt%int_colour  (opt%int_num),
     .                              opt%int_z       (opt%int_num),
     .                              opt%int_a       (opt%int_num),
     .                              opt%int_charge  (opt%int_num),
     .                              opt%int_index   (opt%int_num),
     .                              opt%int_shape   (opt%int_num),
     .                              opt%int_instr   (opt%int_num),
     .                              opt%int_width   (opt%int_num),
     .                              opt%int_database(opt%int_num)
c...            Data source:
                CALL LoadEmissionDatabase(dummy,7,2)

              CASE (3)
c...            Line-of-sight weighted average:
                READ(dummy,*) cdum1,opt%int_type    (opt%int_num),                
     .                              opt%int_colour  (opt%int_num),
     .                              opt%int_z       (opt%int_num),
     .                              opt%int_a       (opt%int_num),
     .                              opt%int_charge  (opt%int_num),
     .                              opt%int_index   (opt%int_num),
     .                              opt%int_average (opt%int_num),
     .                              opt%int_database(opt%int_num)
c...            Data source:
                CALL LoadEmissionDatabase(dummy,7,0)
              CASE DEFAULT
            ENDSELECT

          ELSE
            BACKSPACE 5
            EXIT
          ENDIF
        ENDDO
      ELSE
        BACKSPACE 5
c...    Defaults:
      ENDIF

      WRITE(0,*) 'DONE LOADING INPUT OPTIONS'


      WRITE(0,*) 'TEST:'

c      ALLOCATE(test(20))
c      ALLOCATE(test(1)%test1(20))
c      write(0,*) size( transfer(test,byte) )

      WRITE(0,*) 'INTEGRAL:',opt%int_num
      DO i1 = 1, opt%int_num
        WRITE(0,*) '        :',opt%int_type(i1),opt%int_colour(i1)
        WRITE(0,*) '        :',opt%int_z   (i1),opt%int_charge(i1)
        WRITE(0,*) '        :',opt%int_database(i1),opt%int_line(i1)
      ENDDO
c      STOP 'sdfsd'










      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: LoadDetectorOptions985
c
      SUBROUTINE LoadDetectorOptions985(opt)
      USE mod_out985
      USE mod_out989
      IMPLICIT none

      TYPE(type_test), ALLOCATABLE :: test(:)
      character(len=1), dimension(1) :: byte

      TYPE(type_options985) :: opt     
      TYPE(type_header) :: header
      REAL*8 :: image(1000,1000)
c      REAL*8, TARGET :: image(1000,1000)

      INTEGER i1,ix,iy,fp,idum1,idum2
      CHARACTER cdum1*128,dummy*1024


      TYPE(type_options989) :: opt989

      fp = 5

      READ(fp,'(A256)') dummy
      WRITE(0,*) '>>'//dummy(8:10)//'<<'
      IF   (dummy(8:10).EQ.'Ccd'.OR.dummy(8:10).EQ.'CCD'.OR.
     .      dummy(8:10).EQ.'ccd') THEN
        WRITE(0,*) 'LOADING DETECTOR DATA FROM FILE'
        BACKSPACE 5
        READ(fp,*) cdum1,opt%ccd
        READ(fp,*) cdum1,opt%focallength
        READ(fp,*) cdum1,opt%distortion
        READ(fp,*) cdum1,opt%cen(1:3)
        READ(fp,*) cdum1,opt%roll,opt%tilt,opt%swing
        READ(fp,*) cdum1,opt%width(1:2)
        READ(fp,*) cdum1,opt%angle(1:2)
        READ(fp,*) cdum1,opt%nxbin,opt%nybin
        READ(fp,*) cdum1,opt%sa_nxbin,opt%sa_nybin,opt%sa_opt,
     .                   opt%sa_par1 ,opt%sa_par2
        READ(fp,*) cdum1,opt%fmap

        opt%mask_num = 0
        DO WHILE(.TRUE.)
          READ(fp,'(A1024)') dummy
          IF (dummy(10:13).EQ.'mask') THEN
            opt%mask_num = opt%mask_num + 1
            READ(dummy,*) cdum1,opt%mask_opt(opt%mask_num)

            IF (opt%ccd.EQ.0) CYCLE

            IF     (opt%mask_opt(opt%mask_num).EQ.0) THEN
              opt%mask_num = opt%mask_num - 1               
            ELSEIF (opt%mask_opt(opt%mask_num).EQ.1) THEN
c...          Load image and look for pure/artificial black (pixel=0):
              READ(dummy,*) cdum1,opt%mask_opt(opt%mask_num),
     .                            opt989%imagetype,
     .                            opt989%fimage
              CALL LoadImage(opt989,image,header,1000,1000)                 ! find out image size... as elsewhere... 
              opt%img_opt = 1
              opt%img_nxbin = opt989%nxbin
              opt%img_nybin = opt989%nybin
              DO ix = 1, opt989%nxbin
                DO iy = 1, opt989%nybin
                  opt%img_image(ix,iy) = image(ix,iy)  ! Was using pointer, SUN compiler complained...
                ENDDO
              ENDDO
            ELSEIF (opt%mask_opt(opt%mask_num).EQ.2) THEN
            ELSEIF (opt%mask_opt(opt%mask_num).EQ.3) THEN
            ELSE
              CALL ER('985','Unrecognized mask type',*99)
            ENDIF   
          ELSE
            BACKSPACE 5
            EXIT
          ENDIF
        ENDDO

        READ(fp,*) cdum1 ! views for plotting - yes/no

      ELSE
c...    Defaults:
        opt%ccd = 0
        opt%cen(1) = 2.15D0
        opt%cen(2) = 1.42D0
        opt%cen(3) = 0.0D0
        opt%width(1) = 0.20D-07
        opt%width(2) = 0.20D-07
        opt%angle(1) = 65.0D0
        opt%angle(2) = 65.0D0
c...    Number of rows and columns of pixels:
        opt%nxbin =  11
        opt%nybin =  11           
c...  
        opt%roll = 0.0D0
        opt%tilt = 25.0D0
        opt%swing = 75.0D0
        opt%mask_num = 0

        BACKSPACE 5
      ENDIF


      RETURN
 99   STOP
      END
c
c ======================================================================
c
