!     -*-Mode:f90-*-


MODULE MOD_OUT985

  IMPLICIT NONE

  PRIVATE

!   PUBLIC :: ALLOC_OBJ,DEALLOC_OBJ,ALLOC_PIXEL,DEALLOC_PIXEL
  PUBLIC :: ALLOC_CHORD,DEALLOC_CHORD,ALLOC_SURFACE,ALLOC_VERTEX

  INTEGER    MAX3DSUR  ,MAX3DPTS  ,MAX3DLNK   ,MAX3DVER  ,MAXNPLOT
  PARAMETER (MAX3DSUR=6,MAX3DPTS=4,MAX3DLNK=20,MAX3DVER=8,MAXNPLOT=500)
!  PARAMETER (MAX3DSUR=6,MAX3DPTS=4,MAX3DLNK=8,MAX3DVER=8,MAXNPLOT=500)

  INTEGER    MAX3DSIDE 
  PARAMETER (MAX3DSIDE=6)


  INTEGER, PUBLIC, PARAMETER :: MAXNDET = 20,  &
  &                             MAXNCHORD = 1000, &
  &                             MAXNINT = 10,  &
  &                             ALL_OPTIONS = 1,  &
  &                             DETECTOR_ONLY = 2


  INTEGER, PUBLIC, PARAMETER :: MP_INITIALIZE = 1, MP_INCREASE_SIZE = 2

  INTEGER, PUBLIC, PARAMETER :: VTX_SIZESTEP = 100000
  INTEGER, PUBLIC, SAVE :: nvtx
  INTEGER, PUBLIC, SAVE :: maxvtx
  REAL*8 , PUBLIC, ALLOCATABLE, SAVE :: vtx(:,:)

  INTEGER, PUBLIC, SAVE :: dchord  ! Debug chord index
  LOGICAL, PUBLIC, SAVE :: output  ! Debug status

  TYPE, PUBLIC :: type_surface
    INTEGER :: type                    ! Category: SP_PLANAR_POLYGON=4, SP_LINE_SEGMENT=5
    INTEGER :: isector(2)              ! Toroidal sector index range
    INTEGER :: nvtx
    INTEGER :: ivtx(MAX3DVER)
    REAL*8  :: norm(3)                 ! Surface normal
    REAL    :: count                   ! Adding up rays terminating on the surface -SL,03/04/2019
  ENDTYPE type_surface
  REAL              , PUBLIC, PARAMETER :: SRF_VERSION = 1.0
  INTEGER           , PUBLIC, PARAMETER :: SRF_SIZESTEP = 1000000
  INTEGER           , PUBLIC,              SAVE :: nsrf
  INTEGER           , PUBLIC,              SAVE :: maxsrf
  TYPE(type_surface), PUBLIC, ALLOCATABLE, SAVE :: srf(:)



  TYPE, PUBLIC :: type_intersection
     INTEGER :: obj
     INTEGER :: sur
     INTEGER :: srf
     REAL    :: rdum1
     REAL*8  :: dist
     REAL*8  :: v(3)
  ENDTYPE type_intersection



  TYPE, PUBLIC :: type_3D_object
     INTEGER :: index
     INTEGER :: index_gen          ! Points to generic object group (not yet in use)
     INTEGER :: index_pla          ! Points to plasma object
     INTEGER :: type
     INTEGER :: subtype
     INTEGER :: mode
     INTEGER :: surface
     REAL    :: phi
     INTEGER :: wedge1
     INTEGER :: wedge2
     INTEGER :: colour       ! delete, generic object property
     INTEGER :: orientation  ! delete
     INTEGER :: ik           
     INTEGER :: ir
     INTEGER :: in
     INTEGER :: ivolume
     INTEGER :: sample
     INTEGER :: esurf(MAX3DSIDE)
!    Geometry (original):
     INTEGER :: nsur
     INTEGER :: nver
     INTEGER :: tsur(MAX3DSUR)          ! Type: ...
     INTEGER :: gsur(MAX3DSUR)          ! Geometry type: ...
     INTEGER :: rsur(4,MAX3DSUR)
     INTEGER :: flag(MAX3DSUR)  ! FLAG is temporary...
     INTEGER :: npts(MAX3DSUR)
     INTEGER :: ipts(MAX3DPTS,MAX3DSUR)
     REAL*8  :: v(3,MAX3DVER)  ! Move to a separate list of points rather than storing them with each obj...
!    Geometry (new):
     INTEGER :: nside
     INTEGER :: iside(MAX3DSIDE,2)
!    Reflections:
     INTEGER :: reflec(MAX3DSIDE)         ! delete, not handled by generic object group pointer
     INTEGER :: fudge (MAX3DSIDE)         ! delete...
     REAL    :: factor(MAX3DSIDE)         ! delete...
!      tside,gside, etc. are unchanged from tsur,gsur above, just rename them when finally moving to new geometry 
!    Connection map:
     INTEGER :: nmap(MAX3DSUR)
     INTEGER :: imap(MAX3DLNK,MAX3DSUR)
     INTEGER :: isur(MAX3DLNK,MAX3DSUR)
!    Integration quantities:
     REAL*4  :: quantity(MAXNINT)
     REAL*4  :: path           ! Temporary..?
  ENDTYPE type_3D_object




  TYPE, PUBLIC :: type_view ! Wow, I can really save some space here...
     INTEGER :: index
     INTEGER :: tmp
     INTEGER :: index_spe(MAXNINT)      ! Index of line shape spectrum associated with this pixel
     LOGICAL :: valid
     INTEGER :: xindex
     INTEGER :: yindex
     INTEGER :: type
     INTEGER :: ccd
     INTEGER :: status
     INTEGER :: nxbin         ! Don't store these here, use OPT... clean this up...
     INTEGER :: nybin
     INTEGER :: parent
     INTEGER :: ispec(MAXNINT)
     REAL*8  :: integral(MAXNINT)
     REAL*8  :: average (MAXNINT)
     REAL*8  :: scale
     REAL*8  :: weight
     REAL*8  :: weight_save   ! For tertiary+ reflections
     REAL*8  :: v1(3)
     REAL*8  :: v2(3)
     REAL*8  :: xwidth        ! Don't store these here, use OPT... 
     REAL*8  :: ywidth
     REAL*8  :: xangle
     REAL*8  :: yangle
     REAL*8  :: dxangle
     REAL*8  :: dyangle 
     REAL*8  :: rot(3)
     REAL*8  :: trans(3)
     REAL    :: global_v1(3)
     REAL    :: global_v2(3)
     INTEGER :: otrack
     INTEGER :: ntrack
     INTEGER, POINTER :: tlist(:)    ! INTERGER*2?
     REAL*8 , POINTER :: track(:)    ! REAL*4? 
     INTEGER :: nprofile
     REAL*8 , POINTER :: profile(:,:)  ! REAL*4?   Storing the profile of the sampled emission along the viewing chord
  ENDTYPE type_view

  INTEGER, PUBLIC, PARAMETER :: IT_VWINTER = 1,  &
&                               IT_GBINTER = 2,  &
&                               IT_OBINTER = 3     


  INTEGER, PUBLIC, PARAMETER :: MAXSPECBIN = 100



  INTEGER, PUBLIC, PARAMETER :: MAXIONSTATE = 20  ! Should match the number of possible lines to be integrated

  TYPE, PUBLIC :: type_plasma   
     REAL    :: nD
     REAL    :: nD2
     REAL    :: ne           ! *** IF ADDING THINGS HERE, ALSO UPDATE THE INITIALIZATION IN out985_emission.f ***
     REAL    :: te
     REAL    :: nb
     REAL    :: vb
     REAL    :: tb
     REAL    :: si(MAXNINT)  ! Ionisation state assigned for ion charge state of interest, as regards integration channel
     REAL    :: ni(MAXNINT)
     REAL    :: vi(MAXNINT)  ! 0 - magnitude, ...
     REAL    :: ti(MAXNINT)
     REAL    :: bfield(0:3)  ! B-field: 0 magnitude, 1-3 direction unit vector (x,y,z)
  ENDTYPE type_plasma


  TYPE, PUBLIC :: type_test             ! Checking if this sort of thing is allowed, for IDL transfer...
     INTEGER, ALLOCATABLE :: test1(:)
     INTEGER, ALLOCATABLE :: test2(:)
  ENDTYPE type_test
 


  INTEGER, PARAMETER, PUBLIC :: MAX_OPT_MASK      = 20,  &
  &                             MAX_OPT_MASK_NVTX = 10,  &
  &                             MAX_OPT_RIB       = 50


  TYPE, PUBLIC :: type_element  ! SpaceBall
     REAL      :: version = 1.0
     CHARACTER :: tag*128
     INTEGER   :: n
     INTEGER   :: reflec
     INTEGER   :: colour         
     REAL      :: e(  0:2)
     REAL      :: r(2,0:2)
     REAL      :: h(2,0:2)                       
  ENDTYPE type_element          
  
  TYPE, PUBLIC :: type_options985
     INTEGER   :: load 
     INTEGER   :: ccd
     REAL*8    :: focallength
     REAL*8    :: distortion
     REAL*8    :: cen(3)                     ! xcen,ycen,zcen  ! ...
     REAL*8    :: width(2)                   ! xwidth,ywidth
     REAL*8    :: angle(2)                   ! xangle,yangle
     INTEGER   :: nxbin
     INTEGER   :: nybin                      ! nbin(2) !xbin,ybin
     INTEGER   :: chord_n
     INTEGER   :: chord_opt (  MAXNCHORD)
     REAL*8    :: chord_v1  (3,MAXNCHORD)
     REAL*8    :: chord_v2  (3,MAXNCHORD)
     CHARACTER :: fmap*1024
     INTEGER   :: ndet
     INTEGER   :: det_ccd   (MAXNDET)
     INTEGER   :: det_nxbin (MAXNDET)
     INTEGER   :: det_nybin (MAXNDET)             ! nbin(2) !xbin,ybin
     INTEGER   :: det_istart(MAXNDET)
     INTEGER   :: det_iend  (MAXNDET)
     CHARACTER :: det_fname (MAXNDET)*1024
     INTEGER   :: n
     INTEGER   :: m
     REAL*8    :: roll
     REAL*8    :: tilt
     REAL*8    :: swing
     INTEGER   :: sa_opt
     INTEGER   :: sa_nxbin
     INTEGER   :: sa_nybin
     REAL      :: sa_par1
     REAL      :: sa_par2
     INTEGER   :: obj_nsector
     REAL      :: obj_angle_start
     REAL      :: obj_angle_end
     REAL      :: obj_yrotation
     INTEGER   :: obj_num         ! New object descriptors
     INTEGER   :: obj_type  (40)
     INTEGER   :: obj_option(40)
     INTEGER   :: obj_colour(40)
     INTEGER   :: obj_reflec(40)
     INTEGER   :: obj_fudge (40)
     REAL      :: obj_factor(40)
     INTEGER   :: obj_material(40)
     INTEGER   :: obj_orientation(40)
     INTEGER   :: obj_distort(40)     
     INTEGER   :: obj_n(40,2)
     REAL      :: obj_r(40,2)
     REAL      :: obj_z(40,2)
     REAL      :: obj_psi(40,2)
     REAL*8    :: obj_scale(40)
     REAL      :: obj_yangle(40)
     INTEGER   :: obj_sym(40)
     CHARACTER :: obj_fname(40)*256
     CHARACTER :: obj_slist(40)*256
     CHARACTER :: obj_tag  (40)*256
     INTEGER       :: int_num
     INTEGER       :: int_type(MAXNINT)    ! Put these in a structure?
     INTEGER       :: int_colour(MAXNINT)
     INTEGER       :: int_z(MAXNINT)
     INTEGER       :: int_a(MAXNINT)
     INTEGER       :: int_index(MAXNINT)
     INTEGER       :: int_charge(MAXNINT)
     INTEGER       :: int_shape(MAXNINT)
     REAL          :: int_width(MAXNINT)
     REAL          :: int_instr(MAXNINT)
     INTEGER       :: int_average(MAXNINT)
     INTEGER       :: int_database(MAXNINT)
     CHARACTER(64) :: int_line(MAXNINT) 
     INTEGER       :: int_transition(MAXNINT)         ! PIN
     INTEGER       :: int_component(MAXNINT)
     CHARACTER(80) :: int_adasid(MAXNINT)          ! ADAS
     INTEGER       :: int_adasyr(MAXNINT)
     CHARACTER(3)  :: int_adasex(MAXNINT)
     INTEGER       :: int_isele(MAXNINT)
     INTEGER       :: int_iselr(MAXNINT)
     INTEGER       :: int_iselx(MAXNINT)
     INTEGER       :: int_iseld(MAXNINT)
     REAL          :: int_wlngth(MAXNINT)
     REAL          :: int_wlngthbin(MAXSPECBIN,MAXNINT)
!     INTEGER   :: ob_model        ! OLD object descriptors
!     INTEGER   :: ob_nsector
!     REAL      :: ob_angle_start
!     REAL      :: ob_angle_end
!     REAL      :: ob_yrotation
!     INTEGER   :: ob_stdgrd
!     INTEGER   :: ob_stdgrd_colour
!     INTEGER   :: ob_stdgrd_reflec
!     INTEGER   :: ob_trigrd
!     INTEGER   :: ob_trigrd_colour
!     INTEGER   :: ob_trigrd_reflec
!     INTEGER   :: ob_invgrd
!     INTEGER   :: ob_invgrd_colour
!     INTEGER   :: ob_invgrd_idum1
!     REAL      :: ob_invgrd_xcen
!     REAL      :: ob_invgrd_ycen
!     REAL      :: ob_invgrd_xwidth
!     REAL      :: ob_invgrd_ywidth
!     INTEGER   :: ob_invgrd_nxbin
!     INTEGER   :: ob_invgrd_nybin
!     INTEGER   :: ob_vessel
!     INTEGER   :: ob_wall
!     INTEGER   :: ob_wall_colour
!     INTEGER   :: ob_wall_reflec
!     INTEGER   :: ob_targ
!     INTEGER   :: ob_targ_colour
!     INTEGER   :: ob_targ_reflec
!     INTEGER   :: ob_tube
!     INTEGER   :: ob_tube_colour
!     INTEGER   :: ob_tube_idum1
!     INTEGER   :: ob_line
!     INTEGER   :: ob_line_colour
!     INTEGER   :: ob_line_idum1
!     INTEGER   :: ob_user(10)
!     INTEGER   :: ob_user_colour(10)
!     INTEGER   :: ob_user_reflec(10)
!     INTEGER   :: ob_raw_num
!     INTEGER   :: ob_raw_colour(10)
!     INTEGER   :: ob_raw_reflec(10)
!     INTEGER   :: ob_raw_ind(10,2)
!     INTEGER   :: ob_raw_material(10)
!     INTEGER   :: ob_raw_orientation(10)
!     REAL*8    :: ob_raw_scale(10)
!     CHARACTER :: ob_raw_fname(10)*256
     INTEGER   :: img_opt
     INTEGER   :: img_nxbin
     INTEGER   :: img_nybin
     INTEGER   :: img_nxratio
     INTEGER   :: img_nyratio
     ! Reflection models:
     INTEGER   :: ref_opt
     INTEGER   :: ref_num
     INTEGER   :: ref_model(10)
     INTEGER   :: ref_wlgth(10)
     REAL      :: ref_k(10)
     INTEGER   :: ref_n(10)
     REAL      :: ref_cutoff(10)
     INTEGER   :: ref_ow(10)
     REAL      :: ref_pw(10)
     INTEGER   :: ref_otheta(10)
     REAL      :: ref_dtheta(10)
     INTEGER   :: ref_ophi(10)
     REAL      :: ref_dphi(10)
     
!     REAL*8 :: img_image(1000,1000)
     REAL*8 :: img_image(1100,1100)
!     REAL*8, POINTER :: img_image(1000,1000)
     ! Image masks:
     INTEGER   :: mask_num
     INTEGER   :: mask_opt     (MAX_OPT_MASK)
     INTEGER   :: mask_polarity(MAX_OPT_MASK)
     INTEGER   :: mask_nvtx    (MAX_OPT_MASK)
     INTEGER   :: mask_vtx     (MAX_OPT_MASK,MAX_OPT_MASK_NVTX)
     ! Plots
     INTEGER        :: nplots
     CHARACTER(512) :: plots(MAXNPLOT)
     ! Ribbon grid:
     INTEGER        :: rib_n
     INTEGER        :: rib_option(  MAX_OPT_RIB)
     INTEGER        :: rib_nrad  (  MAX_OPT_RIB)  ! Number of steps between (r1,z1) and (r2,z2)
     REAL           :: rib_r     (2,MAX_OPT_RIB)  ! Line segment defining the origin of the ribbon grid
     REAL           :: rib_z     (2,MAX_OPT_RIB)
     INTEGER        :: rib_nphi  (  MAX_OPT_RIB)
     REAL           :: rib_phi   (2,MAX_OPT_RIB)
     REAL           :: rib_scale (  MAX_OPT_RIB)
     INTEGER        :: rib_trace (  MAX_OPT_RIB)  ! Which field line tracing program to call
     REAL           :: rib_dphi  (  MAX_OPT_RIB)  ! Resolution of the field line tracing step, i.e. how far to go when adjusting the pitch of the trace
     CHARACTER(128) :: rib_tfile (  MAX_OPT_RIB)  ! File name for field line trace data (from CASTEM at the moment)
     INTEGER        :: rib_limit (  MAX_OPT_RIB)  ! Limit the extent of the field line generation
     INTEGER        :: rib_param (4,MAX_OPT_RIB)  
     CHARACTER(128) :: rib_tag   (  MAX_OPT_RIB)

     REAL           :: rib_r1    (  MAX_OPT_RIB)  ! region of interest (OPTION=2)
     REAL           :: rib_r2    (  MAX_OPT_RIB)
     REAL           :: rib_s1    (  MAX_OPT_RIB)
     REAL           :: rib_s2    (  MAX_OPT_RIB)

     INTEGER        :: rib_wipe  (  MAX_OPT_RIB)

  ENDTYPE type_options985


!  TYPE(type_3D_object), PUBLIC, ALLOCATABLE, SAVE :: obj(:)
!  TYPE(type_view)     , PUBLIC, ALLOCATABLE, SAVE :: pixel(:)
  TYPE(type_view)     , PUBLIC, ALLOCATABLE, SAVE :: s_chord(:),p_chord(:)  ! Just for plotting chords... 


!  INTEGER, PUBLIC, SAVE :: NOBJ, NPIXEL, NCHORD
  INTEGER, PUBLIC, SAVE :: nchord,n_pchord

! Object properties:  
  INTEGER,PUBLIC :: OP_INTEGRATION_VOLUME  , OP_EMPTY
  PARAMETER        (OP_INTEGRATION_VOLUME=1, OP_EMPTY=2)

  INTEGER, PUBLIC, PARAMETER :: OP_FLUID_GRID     = 1,  &
  &                             OP_EIRENE_GRID    = 2,  &
  &                             OP_INVERSION_GRID = 3
  

! Surface propertis:
  INTEGER,PUBLIC :: SP_GRID_SURFACE  , SP_GRID_BOUNDARY  , SP_VESSEL_WALL
  PARAMETER        (SP_GRID_SURFACE=1, SP_GRID_BOUNDARY=2, SP_VESSEL_WALL=3)

  INTEGER, PUBLIC, PARAMETER :: SP_PLANAR_POLYGON = 4, SP_LINE_SEGMENT=5

! Surface geometry types:
  INTEGER,PUBLIC :: GT_TD  , GT_TC  
  PARAMETER        (GT_TD=1, GT_TC=2)


  REAL*8,PUBLIC :: D_DEGRAD
  PARAMETER       (D_DEGRAD=1.74532925199D-02)

  INTEGER, PUBLIC :: MAXINTER
  PARAMETER         (MAXINTER=100)


  !...  Utility, wall lists:
  INTEGER, PUBLIC :: nvwlist,ngblist,noblist
  INTEGER, PUBLIC, TARGET, ALLOCATABLE :: vwlist(:,:),gblist(:,:)
  INTEGER, PUBLIC, TARGET :: oblist(40,2)

  INTEGER, PUBLIC :: nvwinter,ngbinter,nobinter
  ! jdemod - these need to be declared allocatable
  TYPE(type_intersection), PUBLIC, ALLOCATABLE :: vwinter(:),gbinter(:),obinter(:)

  INTEGER, PUBLIC :: nelement
  TYPE(type_element), PUBLIC, ALLOCATABLE:: element(:)


CONTAINS


!  SUBROUTINE ALLOC_OBJ(NOBJ)
!    INTEGER, INTENT(IN) :: NOBJ
!    IF (ALLOCATED(OBJ)) THEN  
!       WRITE(0,*) 'DEALLOCATING!'
!       DEALLOCATE (OBJ)
!    ENDIF
!    WRITE(0,*) 'ALLOCATING!',NOBJ
!    ALLOCATE (OBJ(NOBJ))
!    RETURN
!  END SUBROUTINE ALLOC_OBJ

!  SUBROUTINE DEALLOC_OBJ
!    IF (.NOT.ALLOCATED(OBJ)) RETURN
!    DEALLOCATE (OBJ)
!    RETURN
!  END SUBROUTINE DEALLOC_OBJ


!  SUBROUTINE ALLOC_PIXEL(NPIXEL)
!    INTEGER, INTENT(IN) :: NPIXEL
!    IF (ALLOCATED(PIXEL)) RETURN
!    ALLOCATE (PIXEL(NPIXEL))
!    RETURN
!  END SUBROUTINE ALLOC_PIXEL

!  SUBROUTINE DEALLOC_PIXEL
!    IF (.NOT.ALLOCATED(PIXEL)) RETURN
!    DEALLOCATE (PIXEL)
!    RETURN
!  END SUBROUTINE DEALLOC_PIXEL


  SUBROUTINE ALLOC_CHORD(NCHORD)
    INTEGER, INTENT(IN) :: NCHORD
    IF (.NOT.ALLOCATED(s_chord)) ALLOCATE (s_CHORD(NCHORD))
    IF (.NOT.ALLOCATED(p_chord)) ALLOCATE (p_CHORD(NCHORD))
    RETURN
  END SUBROUTINE ALLOC_CHORD

  SUBROUTINE DEALLOC_CHORD  ! *REMOVE/REPLACE*
    IF (ALLOCATED(s_chord)) DEALLOCATE (s_CHORD)
    IF (ALLOCATED(p_chord)) DEALLOCATE (p_CHORD)
    RETURN
  END SUBROUTINE DEALLOC_CHORD


  SUBROUTINE Alloc_Surface(memsize,mode)
    INTEGER :: memsize,mode
    TYPE(type_surface), ALLOCATABLE :: tmpsrf(:)
    INTEGER :: istat,ncall
    SAVE
    IF     (mode.EQ.MP_INITIALIZE) THEN
      IF (ALLOCATED(srf)) THEN
        WRITE(0,*) 'ERROR ALLOC_SURFACE: SRF ARRAY ALREADY ALLOCATED'
        STOP
      ELSE
        ncall = 0
        nsrf = 0
        maxsrf = memsize       
        IF (maxsrf.EQ.-1) maxsrf = SRF_SIZESTEP
        ALLOCATE(srf(maxsrf),STAT=istat)
        IF (istat.NE.0) THEN
          WRITE(0,*) 'ERROR ALLOC_SURFACE: PROBLEM ALLOCATING SRF ARRAY'
          STOP      
        ENDIF
      ENDIF
    ELSEIF (mode.EQ.MP_INCREASE_SIZE) THEN
      IF (.NOT.ALLOCATED(srf)) THEN
        WRITE(0,*) 'ERROR ALLOC_SURFACE: SRF ARRAY NOT ALLOCATED'
        STOP
      ELSE
        ncall = ncall + 1
        WRITE(0,*) 'ALLOC_SURFACE: INCREASING SIZE',nsrf,ncall
        ALLOCATE(tmpsrf(nsrf),STAT=istat)     
        IF (istat.NE.0) THEN
          WRITE(0,*) 'ERROR ALLOC_SURFACE: BAD TMPSRF'
          STOP      
        ENDIF
        tmpsrf(1:nsrf) = srf(1:nsrf)
        DEALLOCATE(srf)
        IF (memsize.EQ.-1) THEN
          maxsrf = maxsrf + SRF_SIZESTEP * ncall
        ELSE
          maxsrf = maxsrf + memsize
        ENDIF
        ALLOCATE(srf(maxsrf))
        IF (istat.NE.0) THEN
          WRITE(0,*) 'ERROR ALLOC_SURFACE: CANNOT REALLOCATE SRF'
          STOP      
        ENDIF
        srf(1:nsrf) = tmpsrf(1:nsrf)
        DEALLOCATE(tmpsrf)
      ENDIF
    ELSE
      WRITE(0,*) 'ERROR ALLOC_SURFACE: UNRECOGNIZED MODE VALUE'
      STOP      
    ENDIF
    srf(1:nsrf)%count = 0.0
    RETURN
  END SUBROUTINE ALLOC_SURFACE


  SUBROUTINE Alloc_Vertex(memsize,mode)
    INTEGER :: memsize,mode
    REAL*8, ALLOCATABLE :: tmpvtx(:,:)
    INTEGER :: istat,i1,ncall
    SAVE
    IF     (mode.EQ.MP_INITIALIZE) THEN
      IF (ALLOCATED(vtx)) THEN
        WRITE(0,*) 'ERROR ALLOC_VERTEX: VTX ARRAY ALREADY ALLOCATED'
        STOP
      ELSE
        ncall = 0
        nvtx = 0
        maxvtx = memsize       
        IF (maxvtx.EQ.-1) maxvtx = VTX_SIZESTEP
        ALLOCATE(vtx(3,maxvtx),STAT=istat)
        IF (istat.NE.0) THEN
          WRITE(0,*) 'ERROR ALLOC_VERTEX: PROBLEM ALLOCATING VTX ARRAY'
          STOP      
        ENDIF
      ENDIF
    ELSEIF (mode.EQ.MP_INCREASE_SIZE) THEN
      IF (.NOT.ALLOCATED(vtx)) THEN
        WRITE(0,*) 'ERROR ALLOC_VERTEX: VTX ARRAY NOT ALLOCATED'
        STOP
      ELSE
        ncall = ncall + 1
        WRITE(0,*) 'ADDVERTEX: INCREASING SIZE',nvtx,ncall
        ALLOCATE(tmpvtx(3,nvtx),STAT=istat)     
        IF (istat.NE.0) THEN
          WRITE(0,*) 'ERROR ALLOC_VERTEX: BAD TMPVTX'
          STOP      
        ENDIF
        DO i1 = 1, 3
          tmpvtx(i1,1:nvtx) = vtx(i1,1:nvtx)
        ENDDO
        DEALLOCATE(vtx)
        IF (memsize.EQ.-1) THEN
          maxvtx = maxvtx + VTX_SIZESTEP * ncall
        ELSE
          maxvtx = maxvtx + memsize
        ENDIF
        ALLOCATE(vtx(3,maxvtx))
        IF (istat.NE.0) THEN
          WRITE(0,*) 'ERROR ALLOC_VERTEX: CANNOT REALLOCATE VTX'
          STOP      
        ENDIF
        DO i1 = 1, 3
          vtx(i1,1:nvtx) = tmpvtx(i1,1:nvtx)
        ENDDO
        DEALLOCATE(tmpvtx)
      ENDIF
    ELSE
      WRITE(0,*) 'ERROR ALLOC_VERTEX: UNRECOGNIZED MODE VALUE'
      STOP      
    ENDIF
    RETURN
  END SUBROUTINE ALLOC_VERTEX



END MODULE MOD_OUT985



MODULE MOD_OUT985_VARIABLES
  USE mod_out985
  IMPLICIT none

  PUBLIC

  INTEGER :: grd_ntorseg  

  TYPE(type_options985) :: opt !, POINTER ::  opt     


  INTEGER :: MAX3D,nobj
  TYPE(type_3D_object), ALLOCATABLE :: obj(:) !, POINTER :: obj(:)


  INTEGER, PUBLIC :: nplasma
  TYPE(type_plasma), ALLOCATABLE :: plasma(:)


  INTEGER, PUBLIC :: nspectrum
  REAL, ALLOCATABLE :: spectrum(:,:)


END MODULE MOD_OUT985_VARIABLES
!
! ======================================================================
!
MODULE MOD_OUT985_PLOTS

  REAL, PUBLIC :: map1x,map2x,map1y,map2y

END MODULE MOD_OUT985_PLOTS
!
! ======================================================================
!
MODULE mod_out985_ribbon
  USE mod_out985

  PUBLIC

  INTEGER :: trace_n  ,  &
             maxntrace,  &
             ntrace   ,  &
             wipe_n   ,  &
             wipe_list(100)

  INTEGER, ALLOCATABLE :: trace_i(:,:)
  REAL*8 , ALLOCATABLE :: trace  (:,:)

  REAL*8 crop_r1,crop_r2,crop_s1,crop_s2  ! active region of interest

  CONTAINS
! ----------------------------------------------------------------------
  SUBROUTINE AllocTrace(memsize,mode)
    IMPLICIT none

    INTEGER, INTENT(IN) :: memsize,mode

    INTEGER, PARAMETER :: TRACE_SIZESTEP = 100000

    REAL*8, ALLOCATABLE :: tmptrace(:,:)
    INTEGER :: istat,i1

    SELECTCASE (mode)
      ! ----------------------------------------------------------------
      CASE (MP_INITIALIZE) 
        IF (ALLOCATED(trace)) THEN
          WRITE(0,*) 'ERROR AllocTrace: TRACE array already allocated'
          STOP
        ENDIF
        ntrace = 0
        maxntrace = memsize       
        IF (maxntrace.EQ.-1) maxntrace = TRACE_SIZESTEP
        ALLOCATE(trace(5,maxntrace),STAT=istat)
        IF (istat.NE.0) THEN
          WRITE(0,*) 'ERROR AllocTrace: Problem allocating array'
          STOP      
        ENDIF
        trace = 0.0D0
      ! ----------------------------------------------------------------
      CASE (MP_INCREASE_SIZE)
        IF (.NOT.ALLOCATED(trace)) THEN
          WRITE(0,*) 'ERROR AllocTrace: Cannot resize, array not allocated'
          STOP
        ENDIF
        ALLOCATE(tmptrace(5,ntrace),STAT=istat)     
        IF (istat.NE.0) THEN
          WRITE(0,*) 'ERROR AllocTrace: Cannot allocate TMPTRACE'
          STOP      
        ENDIF
        tmptrace = 0.0D0
        DO i1 = 1, 5
          tmptrace(i1,1:ntrace) = trace(i1,1:ntrace)
        ENDDO
        DEALLOCATE(trace)
        IF (memsize.EQ.-1) THEN
          maxntrace = maxntrace + TRACE_SIZESTEP
        ELSE
          maxntrace = maxntrace + memsize
        ENDIF
        ALLOCATE(trace(5,maxntrace))
        IF (istat.NE.0) THEN
          WRITE(0,*) 'ERROR AllocTrace: Cannot reallocate TRACE'
          STOP      
        ENDIF
        DO i1 = 1, 5
          trace(i1,1:ntrace) = tmptrace(i1,1:ntrace)
        ENDDO
        DEALLOCATE(tmptrace)
      ! ----------------------------------------------------------------
      CASE DEFAULT
        WRITE(0,*) 'ERROR AllocTrace: Unrecognised mode'
        STOP      
      ! ----------------------------------------------------------------
    ENDSELECT
    RETURN
  END SUBROUTINE AllocTrace
!
! ----------------------------------------------------------------------
!
END MODULE mod_out985_ribbon
!
! ======================================================================
!
MODULE MOD_OUT985_CLEAN

  PUBLIC :: ALLOC_RESET

  CONTAINS

  SUBROUTINE DEALLOC_ALL
    USE mod_out985
    USE mod_out985_variables
    USE mod_out985_ribbon
    nobj = 0
    nsrf = 0
    nvtx = 0
    ntrace = 0
    IF (ALLOCATED(obj  )) DEALLOCATE (obj  )
    IF (ALLOCATED(srf  )) DEALLOCATE (srf  )
    IF (ALLOCATED(vtx  )) DEALLOCATE (vtx  )
    IF (ALLOCATED(trace)) DEALLOCATE (trace)
    RETURN
  END SUBROUTINE DEALLOC_ALL

END MODULE MOD_OUT985_CLEAN


