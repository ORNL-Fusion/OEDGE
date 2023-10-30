c     -*Fortran*-
c
c ======================================================================
c
c  subroutine: WriteEireneFiles_06
c
c
c
c
c
      SUBROUTINE WriteEireneFiles_06(iitersol)
      USE mod_eirene06
      USE mod_eirene06_locals
      USE mod_sol28_global
      USE mod_options
      USE mod_geometry
      USE mod_filament
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_slcom
      use debug_options
      IMPLICIT none

c     include 'params'
c     include 'cgeom'
c     include 'comtor'
c     INCLUDE 'slcom'

      INTEGER, INTENT(IN) :: iitersol

      INTEGER   ik,ir,in1,in2,i1,id,ik1,status,i,ivoid,istart,iend,
     .          tetrahedron_source,itet, nvoid_store,
     .          fp,idum1,idum2
      LOGICAL   saved_triangles,output
      CHARACTER fname*1024

      REAL t  ! *** TEMP ***

      DATA nvoid_store /-1/
      DATA t /0.0/
      DATA saved_triangles /.FALSE./
      SAVE

      call pr_trace('WriteEireneFiles_06','START')
      
      IF (opt_fil%opt.NE.0) THEN
        IF (citersol.GT.0) THEN
          IF (t.EQ.0.0) THEN
            t = 1.0E-07
c            t = t - opt_fil%time_step + 1.0E-09  ! The last one is so that t=0.0 is registered again  
            CALL DefineFilaments 
            DO i1 = 1, 1000
              CALL SetupFilaments(DBLE(t))
              IF ((t-opt_fil%start_time).GT.-1.0E-06) EXIT
              t = t + opt_fil%time_step
            ENDDO
          ELSE
            t = t + opt_fil%time_step
            CALL SetupFilaments(DBLE(t))
          ENDIF
          WRITE(0,*) 
          WRITE(0,*) '------------------------------------------'
          WRITE(0,*) ' FILAMENT SETUP:',iitersol,NINT(t/1.0E-06),
     .                                  nfilament
          WRITE(0,*) '------------------------------------------'
          WRITE(0,*) 
        ELSE
          IF (t.EQ.0.0) THEN
            CALL DefineFilaments 
            DO i1 = 1, 1000
              CALL SetupFilaments(DBLE(t))
              IF ((t-opt_fil%start_time).GT.-1.0E-06) EXIT
              t = t + opt_fil%time_step
            ENDDO
            WRITE(0,*) 'NON-ITERATIVE TIME STAMP=',t/1.0E-06
          ENDIF
        ENDIF
      ENDIF

      helium = .FALSE.

      output = .FALSE.
      !output = .TRUE.

      opt%pin_data = .TRUE.
      opt_iteration(1:nopt)%pin_data = .TRUE.

      tetrahedrons = .FALSE.
      time_dependent = .FALSE.

      SELECTCASE (eirgeom)
        CASE (-1)
          CALL ER('WriteEireneFiles_06','Bad EIRGEOM',*99)
        CASE ( 3)
          tetrahedrons = .TRUE.
      ENDSELECT        

      IF (opt_eir%ntime.NE.0) THEN
        time_dependent = .TRUE.
        IF (iitersol.EQ.1) THEN
          time0 = opt_eir%time0
        ELSE
          time0 = time0 + opt_eir%dtimv
        ENDIF
      ENDIF

c      IF (tetrahedrons) THEN
c...    All a bit cavalier at the moment because the object parameters
c       from the input file are not currently saved, so can't
c       do a proper check on whether this grid is the right one...
c        CALL LoadObjects('tetrahedrons.raw',status)
c        IF (status.NE.-1) savedtriangles = .TRUE.
c      ENDIF

      CALL OutputData(86,'Before calling EIRENE')

c...  Some variables used in the routines to write the '06 files, which do not have
c     access to the standard DIVIMP common blocks:

      fluid_code = 'DIVIMP'

      time  = eirtime                           ! EIRENE run time (s)
      niter = eirniter                          ! Number of EIRENE self-iterations
      dtimv = opt_eir%dtimv                     ! Time-dependent runs

      ttemp = ctargt                            ! Temperature (K) of target surfaces
      wtemp = cwallt                            ! Temperature (K) of wall surfaces

c...  Surface material, for use with the local reflection models
c     specified in block 3 in the EIRENE input file:
      IF (eirmat1.EQ.1) tmater = 9642.0         ! Target material
      IF (eirmat1.EQ.2) tmater = 1206.0
      IF (eirmat1.EQ.3) tmater = 18474.0
      IF (eirmat1.EQ.4) tmater = 904.0  
      IF (eirmat1.EQ.5) tmater = 5626.0
      IF (eirmat2.EQ.1) wmater = 9642.0         ! Wall
      IF (eirmat2.EQ.2) wmater = 1206.0
      IF (eirmat2.EQ.3) wmater = 18474.0
      IF (eirmat2.EQ.4) wmater = 904.0
      IF (eirmat2.EQ.5) wmater = 5626.0 

      opacity    = eiropacity
      photons    = eirphoton
      trim_data  = eirtrim
      ermin      = eirermin
      bgk        = eirbgk
      ntorseg    = eirntorseg
      torfrac    = eirtorfrac
      alloc      = eiralloc
      whipe      = opt_eir%whipe
      fluid_grid = opt_eir%fluid_grid
      beam       = 0

      i1trc = opt_eir%i1trc 
      i2trc = opt_eir%i2trc 

      eirfp = 88

      eir_pass = 0

      IF (eirtrim.NE.1) CALL ER('WriteEireneFiles','TRIM option '//
     .                          'not set',*99)

      torus1 = 0.0
      torus2 = eirtorfrac * 360.0 * 0.99999  !190.0 ! 180.0

c...  A couple of checks:
      IF (eirniter.EQ.0.AND.(photons.GT.1.OR.bgk.EQ.3))
     .  CALL ER('WriteEireneFiles','Bad EIRNITER',*99)
      IF (eiropacity.NE.0.AND.eirphoton.NE.0)
     .  CALL WN('WriteEireneFiles','Bad EIRNITER')
c     .  WRITE(0,*) 'WARNING: EIROPACITY.NE.0 AND EIRPHOTON.NE.0'

c...  Main calls for generating Eirene geometry objects (triangles 
c     and tetrahedrons at the moment):
      IF (.TRUE.) THEN

        CALL ALLOC_CELL(MAXNKS*nrs,nrs*2)

c        IF (.true.) THEN
c        IF (.false.) THEN
        IF (opt_eir%f_eirene_load.NE.0.OR.saved_triangles) THEN
c        IF (saved_triangles.AND..NOT.tetrahedrons) THEN

          IF (opt_eir%f_eirene_load.NE.0.AND.
     .        .NOT.saved_triangles) CALL CollectEireneFiles

          IF (output) WRITE(0,*) 'building triangles'
          CALL ALLOC_VERTEX  (50000)   ! NEED TO ADD ACTIVE BOUNDS CHECKING FOR ALL THESE!
          CALL ALLOC_SURFACE (3000)
          CALL ALLOC_TRIANGLE(100000)
          ntri = 0
          ntry = 0
          nvtx = 0


          CALL DefineEireneSurfaces_06

          CALL ProcessFluidGrid_06

          CALL LoadTriangles_06

          CALL AssignPlasmaQuantities_06
          CALL SetupEireneStrata

          IF (tetrahedrons) THEN          

            IF (opt_fil%opt.NE.0) THEN
              WRITE(0,*) 'CODE NOT READY FOR FILAMENTS'
              STOP
            ENDIF

            CALL LoadObjects('tetrahedrons.raw',status)
            write(0,*) 'loaded objects',nobj,nsrf,nvtx

c            CALL BuildConnectionMap(1,nobj)
c            CALL RemoveDuplicateVertices

            IF (ALLOCATED(plasma)) DEALLOCATE(plasma) ! not necessary?
            IF (ALLOCATED(bfield)) DEALLOCATE(bfield) ! not necessary?
            IF (ALLOCATED(vtxmap)) DEALLOCATE(vtxmap) ! necessary -- not sure why...
            fp = 99
            OPEN(UNIT=fp,FILE='tetrahedrons_sup.raw',
     .           ACCESS='SEQUENTIAL',
     .           FORM='UNFORMATTED',STATUS='OLD',ERR=98)            
            READ(fp,ERR=98) idum1,nvtxmap,idum2

            write(0,*) 'ntry',idum2
            ALLOCATE(plasma(20,idum2))
            ALLOCATE(bfield(4 ,idum2))
            ALLOCATE(vtxmap(   nvtx))
            READ(fp,ERR=98) vtxmap
            READ(fp,ERR=98) plasma
            READ(fp,ERR=98) bfield
            CLOSE (fp)      

          ENDIF

          saved_triangles = .TRUE.

        ELSE
          IF (output) WRITE(0,*) 'building triangles'
          CALL ALLOC_VERTEX  (50000)   ! NEED TO ADD ACTIVE BOUNDS CHECKING FOR ALL THESE!
          CALL ALLOC_SURFACE (3000)
          CALL ALLOC_TRIANGLE(100000)

          CALL DefineEireneSurfaces_06

c         (The calls below can be reordered, eventually, to make things a little more efficient.)

c...      Loads cell geometry and plasma quantities from fluid grid:
          CALL ProcessFluidGrid_06

c...      Generates Eirene triangle mesh from fluid grid information:
          CALL BuildFluidGridTriangles_06

c...      Builds connection map, amongst other things:
          CALL ProcessTriangles_06(0)

c...      Sets triangle volume quantities from data stored in ProcessFluidGrid:
          CALL AssignPlasmaQuantities_06

c...      Define Eirene particle sources:
          CALL SetupEireneStrata

c...      Fills the voids outside the fluid grid with triangles (by calling
c         the external program TRIANGLE):
          IF (opt_eir%nvoid.GT.0) THEN
            IF (nvoid_store.EQ.-1) THEN  ! hack, since nvoid is increased during processing, which is breaks iteration
              nvoid_store = opt_eir%nvoid
            ELSE
              opt_eir%nvoid = nvoid_store 
            ENDIF
            CALL SetupVoidProcessing(opt_eir)
            DO ivoid = 1, opt_eir%nvoid
              IF (opt_eir%void_version.EQ.1.0) THEN
                CALL ProcessVoid_v1_0(opt_eir%void_zone(ivoid),opt_eir)
              ELSE
                CALL ProcessVoid     (opt_eir%void_zone(ivoid),opt_eir)
              ENDIF
            ENDDO
          ELSE
            CALL WritePolyFile_06(eirntri,MAXNRS,eirtri)
          ENDIF

c...      Re-builds connection map, taking into account the new triangles:
          CALL ProcessTriangles_06(0)

c...      Dumps trianlges to a binary file, for use with LoadTriangles:
          CALL SaveTriangles_06
          saved_triangles = .TRUE.

c...      Tetrahedrons:
          IF (tetrahedrons) THEN
            tetrahedron_source = 1
            DO itet = 1, opt_eir%tet_n
              IF (opt_eir%tet_type(itet).GE.2.0.AND.
     .            opt_eir%tet_type(itet).LE.4.0) THEN
                tetrahedron_source = 2
                EXIT
              ENDIF
            ENDDO

            IF (tetrahedron_source.EQ.1) THEN
              CALL ProcessTetrahedrons_06
            ELSE
              CALL AssembleTetrahedrons
            ENDIF

            IF (opt_fil%opt.NE.0.AND.citersol.GT.0) THEN 
              WRITE(fname,'(A,I3.3,A)') '.',iitersol,'.raw'
              WRITE(0,*) 'FILENAME=','tetrahedrons'//TRIM(fname)
              CALL SaveGeometryData('tetrahedrons'//TRIM(fname))
              CALL SaveFilamentData('filaments'//TRIM(fname))
            ELSE
              CALL SaveGeometryData('tetrahedrons.raw')
              CALL SaveFilamentData('filaments.raw')
            ENDIF
c            CALL DumpGrid('BUILDING TETRAHEDRONS')          
          ELSE
            IF (fluid_grid.EQ.2.OR.fluid_grid.EQ.3) 
     .        CALL AssignTriangleDumpSurface  ! fluid triangles
          ENDIF

c          WRITE(0,*) 'DONE'
        ENDIF

        IF (.TRUE.) CALL LocalGridRefinement

c...    Writes the .points, .sides, .map and .plasma files to be loaded
c       by EIRENE:

c        write(0,*) 'Writing EIRENE objects',tetrahedrons

        IF (tetrahedrons) THEN
          CALL WriteEireneObjects
        ELSE        
          CALL WriteEireneTriangles
        ENDIF

        IF (eirdata.EQ.1) CALL WriteEireneInputFile_06
      ENDIF

c...  Fluid grid data no longer required:
      CALL DEALLOC_CELL

c      CALL DumpGrid('Quitting early')  
c      STOP 'QUITTING EARLY...'

c...  Need this to generate celldata.dat for triangle grid, which
c     is required for finding the B-field everywhere in the vessel:
c      CALL ProcessTriangles

c      IF (tetrahedrons) CALL DumpGrid('BUILDING TETRAHEDRONS')  



      IF (opt_eir%gauge_n.GT.0) CALL eirWritePressureGaugeFile  



      RETURN
 98   WRITE(0,*) '_sup file not found'
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE CollectEireneFiles
      USE mod_sol28
      USE mod_options
      IMPLICIT none

      INTEGER        status
      CHARACTER*1024 fname,command


      write(0,*) 'collecting eirene files'

c...  
      fname = TRIM(opt_eir%f_eirene_15)//'.raw.tri'

c...  
      command = 'cp '//TRIM(opt_eir%f_eirene_dir)//TRIM(fname)//'.gz .'
      write(0,*) '  command: '//TRIM(command)
      CALL CIssue(TRIM(command),status)
      IF (status.NE.0) 
     .  CALL ER('...','Unable to copy file',*99)

      CALL UnzipFile(fname//'.gz')

c...  
      command = 'mv '//TRIM(fname)//' triangles.raw'
      write(0,*) '  command: '//TRIM(command)
      CALL CIssue(TRIM(command),status)
      IF (status.NE.0) 
     .  CALL ER('...','Unable to rename file',*99)


c...  
      fname = TRIM(opt_eir%f_eirene_15)//'.raw.tet'

c...  
      command = 'cp '//TRIM(opt_eir%f_eirene_dir)//TRIM(fname)//'.gz .'
      write(0,*) '  command: '//TRIM(command)
      CALL CIssue(TRIM(command),status)
      IF (status.NE.0) 
     .  CALL ER('...','Unable to copy file',*99)

      CALL UnzipFile(fname//'.gz')

c...  
      command = 'mv '//TRIM(fname)//' tetrahedrons.raw'
      write(0,*) '  command: '//TRIM(command)
      CALL CIssue(TRIM(command),status)
      IF (status.NE.0) 
     .  CALL ER('...','Unable to rename file',*99)



c...  
      fname = TRIM(opt_eir%f_eirene_15)//'.raw.tet_sup'

c...  
      command = 'cp '//TRIM(opt_eir%f_eirene_dir)//TRIM(fname)//'.gz .'
      write(0,*) '  command: '//TRIM(command)
      CALL CIssue(TRIM(command),status)
      IF (status.NE.0) 
     .  CALL ER('...','Unable to copy file',*99)

      CALL UnzipFile(fname//'.gz')
c...  
      command = 'mv '//TRIM(fname)//' tetrahedrons_sup.raw'
      write(0,*) '  command: '//TRIM(command)
      CALL CIssue(TRIM(command),status)
      IF (status.NE.0) 
     .  CALL ER('...','Unable to rename file',*99)





      RETURN
 99   WRITE(0,*) '  FILE NAME = ',TRIM(fname)
      WRITE(0,*) '  COMMAND   = ',TRIM(command)
      WRITE(0,*) '  ERROR     = ',status
      STOP
      END
c
c ======================================================================
c
c subroutine: eirWritePressureGaugeFile
c
      SUBROUTINE eirWritePressureGaugeFile
      USE mod_options
      USE mod_eirene_history
      use mod_params
      IMPLICIT none

c     include 'params'

      INTEGER       :: fp,n,i1,i2,count
      REAL          :: x,y,z,phi      
      CHARACTER*128 :: tag

      LOGICAL, SAVE :: initialized = .FALSE.



c...  Only need to call this routine once:
      IF (initialized) RETURN

      initialized = .TRUE.



      fp = 99
      OPEN(UNIT=fp,FILE='eirene.gauges',ACCESS='SEQUENTIAL',
     .     STATUS='NEW',ERR=97)

      n = 0
      DO i1 = 1, opt_eir%gauge_n
        n = n + opt_eir%gauge_dupe_n(i1) + 1
      ENDDO
      WRITE(fp,*) n

      count = 0
      phi   = -999.0
      DO i1 = 1, opt_eir%gauge_n
        DO i2 = 1, opt_eir%gauge_dupe_n(i1) + 1


          IF (opt_eir%gauge_dupe_n(i1).GT.0) THEN
            WRITE(tag,'(A,I8)') TRIM(opt_eir%gauge_tag(i1)),i2
          ELSE
            tag = TRIM(opt_eir%gauge_tag(i1))
          ENDIF

          IF (i2.EQ.1) THEN

            phi = opt_eir%gauge_phi(i1)

          ELSE

            SELECTCASE (TRIM(opt_eir%gauge_dupe_dir(i1)))
c             ----------------------------------------------------------
              CASE ('tor')  
                phi = opt_eir%gauge_phi      (i1) - REAL(i2-1) * 
     .                opt_eir%gauge_dupe_step(i1)
c             ----------------------------------------------------------
              CASE DEFAULT
                CALL ER('eirWritePressureGaugeFile','Unknown '//
     .                  'duplication direction',*99)
c             ----------------------------------------------------------
            ENDSELECT

          ENDIF

c...      Rotate the vector:
          phi = phi * DEGRAD 
          x = COS(phi) * opt_eir%gauge_x(i1) - SIN(phi) * 0.0
          y =            opt_eir%gauge_y(i1)
          z = SIN(phi) * opt_eir%gauge_x(i1) + COS(phi) * 0.0

 
          WRITE(fp,'(4F12.4,5X,A)') 
     .      x*100.0,y*100.0,z*100.0,opt_eir%gauge_radius(i1)*100.0,tag

          count = count + 1
          opt_eir%gauge_ind(1:2,count) = (/ i1, i2 /)
          opt_eir%gauge_pos(1:4,count) = (/ x, y, z, phi /)

          write(0,*) 'gauge',i1,phi/DEGRAD,x,y,z,tag
 
        ENDDO
      ENDDO

      CLOSE(fp)


      RETURN
 97   CALL ER('eirWritePressureGaugeFile','Could not open file',*99)
 99   STOP
      END
c
c ======================================================================
c
c subroutine: DefineEireneSurfaces
c
      SUBROUTINE DefineEireneSurfaces_06
      USE mod_interface
      USE mod_eirene06_parameters
      USE mod_eirene06
      USE mod_sol28_io
      USE mod_sol28_global
      use mod_params
      use mod_comtor
      use mod_cgeom
      use mod_slcom
      IMPLICIT none
c     INCLUDE 'params' 
c     INCLUDE 'comtor'
c     INCLUDE 'cgeom'
c     INCLUDE 'slcom'

      INTEGER NewEireneSurface_06
      LOGICAL CheckIndex,osmGetLine,osmCheckTag

      INTEGER i1,i2,ik,ik1,ik2,side,sur1,sur2,sur3,iliin,ntmp,
     .        type,index1,index2,ilside,ilswch,region,code,is,nboundary,
     .        tmp_ilspt,fp,n,surface_type,hole_index
      LOGICAL assigned(2),first_message,active
      REAL    x1,x2,xcen,y1,y2,z1,z2,ycen,angle,dangle,rad,ewall,
     .        mater,recycf,recyct,ilspt,isrs,recycs,recycc,
     .        transp1,transp2,file_version
      CHARACTER buffer*1024,fname*512,ftag*128

      first_message = .TRUE.

      nsurface = 0
      default_surface = 0
      core_boundary   = 0

      opt_eir%ilspt = 0

c     ----------------------------------------------------------------------
c...  Setup non-default standard surfaces related to the magnetic grid:
c     ----------------------------------------------------------------------

c...  Poloidal surfaces for Eirene strata:
      IF (.TRUE.) THEN
c...    Determine the number of target strata that are defined:
c        assigned = .FALSE.
c        DO is = 1, opt_eir%nstrata
c          IF (NINT(opt_eir%type(is)) .EQ.1.AND.
c     .             opt_eir%target(is).EQ.IKLO) assigned(IKLO) =.TRUE.
c          IF (NINT(opt_eir%type(is)) .EQ.1.AND.
c     .             opt_eir%target(is).EQ.IKHI) assigned(IKHI) =.TRUE.
c        ENDDO

        nboundary = 0

        assigned = .TRUE.
c       --------------------------------------------------------------------
c...    Define the default surfaces for the target strata, since target
c       strata are not specifically assigned:
        DO i1 = IKLO, IKHI 
          IF (assigned(i1)) CYCLE 
          STOP 'TURING THIS OBSOLETE STRATA CODE OFF!'
          nsurface = NewEireneSurface_06(NON_DEFAULT_STANDARD)
          surface(nsurface)%subtype  = STRATUM
          surface(nsurface)%surtxt   = '* default target (DIVIMP)'
          nboundary = nboundary - 1
          IF (cgridopt.EQ.LINEAR_GRID.OR.cgridopt.EQ.RIBBON_GRID) THEN
            surface(nsurface)%index(1) = irsep               ! Ring index start location of surface
            surface(nsurface)%index(2) = nrs-1               ! Ring index end
            surface(nsurface)%index(3) = i1                  ! Target (IKLO=inner, IKHI=outer)
            surface(nsurface)%index(5) = 0
            surface(nsurface)%index(6) = -nboundary
          ELSE
            surface(nsurface)%index(1) = irsep               ! Ring index start location of surface
            surface(nsurface)%index(2) = nrs                 ! Ring index end
            surface(nsurface)%index(3) = i1                  ! Target (IKLO=inner, IKHI=outer)
            surface(nsurface)%index(5) = 0
            surface(nsurface)%index(6) = -nboundary
          ENDIF
          surface(nsurface)%reflect = LOCAL                  ! Set surface reflection model to LOCAL
          surface(nsurface)%iliin  = 1
          surface(nsurface)%ilside = 0
          surface(nsurface)%ilswch = 0
          surface(nsurface)%iltor  = 0  
          surface(nsurface)%ilcell = 0
          surface(nsurface)%ilcol  = 4
          surface(nsurface)%material = tmater                ! Set surface material
          surface(nsurface)%ewall = -ttemp * 1.38E-23 / ECH  ! Set temperature
          surface(nsurface)%ilspt = 0 ! opt_eir%ilspt(is)
        ENDDO
c...    Loop over the user specified strata and assemble the 
c       corresponding target surfaces:
        n = opt_eir%nstrata
        CALL inOpenInterface('idl.eirene_strata',ITF_WRITE)
        CALL inPutData(opt_eir%type      (  1:n),'TYPE'  ,'N/A')
        CALL inPutData(opt_eir%target    (  1:n),'TARGET','N/A')
        CALL inPutData(opt_eir%range_tube(1,1:n),'RANGE1','N/A')
        CALL inPutData(opt_eir%range_tube(2,1:n),'RANGE2','N/A')
        CALL inCloseInterface
        DO is = 1, opt_eir%nstrata 
          IF (NINT(opt_eir%type(is)).NE.1) CYCLE

          nboundary = nboundary - 1
          nsurface = NewEireneSurface_06(NON_DEFAULT_STANDARD)
          surface(nsurface)%subtype  = STRATUM
          surface(nsurface)%surtxt   ='* user specified target (DIVIMP)'
          surface(nsurface)%index(1:2) = opt_eir%range_tube(1:2,is) ! irsep                  ! Ring index start location of surface
c          surface(nsurface)%index(2) = nrs                    ! Ring index end
          surface(nsurface)%index(3) = opt_eir%target(is)  ! Target (IKLO=inner, IKHI=outer)
          surface(nsurface)%index(5) = 0
          surface(nsurface)%index(6) = nboundary
          surface(nsurface)%reflect = LOCAL                   ! Set surface reflection model to LOCAL
          surface(nsurface)%iliin  = 1
          surface(nsurface)%ilside = 0
          surface(nsurface)%ilswch = 0
          surface(nsurface)%iltor  = 0  
          surface(nsurface)%ilcell = 0
          surface(nsurface)%ilcol  = 4
          surface(nsurface)%material = tmater                ! Set surface material
          surface(nsurface)%ewall = -ttemp * 1.38E-23 / ECH  ! Set temperature
          surface(nsurface)%ilspt = 0 ! opt_eir%ilspt(is)
        ENDDO
      ELSE
        STOP 'OBSOLETE STRATUM SETUP CODE'
      ENDIF

c     ----------------------------------------------------------------------
c...  Core boundary surface:
      nsurface = NewEireneSurface_06(NON_DEFAULT_STANDARD)
      core_boundary = nsurface
      surface(nsurface)%subtype  = MAGNETIC_GRID_BOUNDARY
      nboundary = 1
      surface(nsurface)%index(1) = 1       ! Cell index start location of surface
      surface(nsurface)%index(2) = MAXNKS  ! Cell index end   
      surface(nsurface)%index(3) = 2       ! Ring in the magnetic grid 
      surface(nsurface)%index(4) = 14      ! Radial surface of cell (1-4 cell side here)
      surface(nsurface)%index(5) = 0
      surface(nsurface)%index(6) = nboundary
      surface(nsurface)%ilcol  = 3         ! Colour of surface in EIRENE plots
      IF (cgridopt.EQ.LINEAR_GRID) THEN
        surface(nsurface)%surtxt   = '* core, reflecting (DIVIMP)'
        surface(nsurface)%iliin  = 1         
      ELSE
        surface(nsurface)%surtxt   = '* core, absorbing (DIVIMP)'
        surface(nsurface)%iliin  = 2         ! Reflection type (ILIIN=2 is 100% absorbing)
      ENDIF

c     ----------------------------------------------------------------------
c...  SOL radial boundary surfaces:
      ik1 = 0
      ik2 = 0
      DO ik = 1, nks(irwall)
c      DO ik = 1, nks(irwall)-1
c...
        IF (ik1.EQ.0) ik1 = ik
        IF (ik.LT.nks(irwall).AND.
     .      ((ikins(ik,irwall).NE.ikins(ik+1,irwall)-1).OR.
     .       (irins(ik,irwall).NE.irins(ik+1,irwall)  ))) ik2 = ik
        IF (ik.EQ.nks(irwall).AND.ik2.EQ.0) ik2 = nks(irwall)
c        IF ((ikins(ik,irwall).NE.ikins(MIN(nks(irwall),ik+1),irwall)-1).OR.
c     .      (irins(ik,irwall).NE.irins(ik+1,irwall)  )) ik2 = ik
c        IF (ik.EQ.nks(irwall)-1.AND.ik2.EQ.0) ik2 = nks(irwall)

        IF (ik2.NE.0) THEN
          IF (irouts(ikins(ik1,irwall),irins(ik1,irwall)).EQ.
     .        irwall) THEN
            side = 23
          ELSE
            side = 14
          ENDIF
          nsurface = NewEireneSurface_06(NON_DEFAULT_STANDARD)
          nboundary = nboundary + 1
          surface(nsurface)%subtype  = MAGNETIC_GRID_BOUNDARY
          surface(nsurface)%surtxt   = '* radial SOL (DIVIMP)'
          surface(nsurface)%index(1) = ikins(ik1,irwall)  ! Cell index start location of surface       
          surface(nsurface)%index(2) = ikins(ik2,irwall)  ! Cell index end                             
          surface(nsurface)%index(3) = irins(ik1,irwall)  ! Ring in the magnetic grid                  
          surface(nsurface)%index(4) = side               ! Radial surface of cell
          surface(nsurface)%index(5) = 0                                                                
          surface(nsurface)%index(6) = nboundary                                                        
          surface(nsurface)%iliin  = -1                   ! Reflection type (ILIIN=-1 is 100% transparent)
          surface(nsurface)%ilcol  = 4                    ! Colour of surface in EIRENE plots
          ik1 = 0
          ik2 = 0
        ENDIF
      ENDDO

c     ----------------------------------------------------------------------
c...  PFZ radial boundary:  
      IF (cgridopt.EQ.LINEAR_GRID.OR.irtrap.GT.nrs.OR.
     .    cgridopt.EQ.RIBBON_GRID) THEN
      ELSE
        nboundary = nboundary + 1
        nsurface = NewEireneSurface_06(NON_DEFAULT_STANDARD)
        surface(nsurface)%subtype  = MAGNETIC_GRID_BOUNDARY
        surface(nsurface)%surtxt   = '* radial PFR (DIVIMP)'
        surface(nsurface)%index(1) = 1                         ! Cell index start location of surface       
        surface(nsurface)%index(2) = MAXNKS                    ! Cell index end                             
        surface(nsurface)%index(3) = irtrap + 1                ! Ring in the magnetic grid                  
        surface(nsurface)%index(4) = 14                        ! Radial surface of cell (1-4 cell side here)
        surface(nsurface)%index(5) = 0                                                                      
        surface(nsurface)%index(6) = nboundary                                                        
        surface(nsurface)%iliin  = -1                          ! Reflection type (ILIIN=-1 is 100% transparent)
        surface(nsurface)%ilcol  = 4                           ! Colour of surface in EIRENE plots
      ENDIF

c     ------------------------------------------------------------------
c     Add wall surfaces to the list (not for inclusion in the EIRENE 
c     input file):

c...  Load vessel wall segments from the default DIVIMP vessel 
c     wall array (WALLPTS), clockwise wall is assumed:
      DO i1 = wallpts, 1, -1
        IF (i1.EQ.wallpts) THEN
          i2 = 1
        ELSE
          i2 = i1 + 1
        ENDIF
        IF (wallpt(i1,18).EQ.0.0) THEN
          nsurface = NewEireneSurface_06(VESSEL_WALL)
          surface(nsurface)%index(1) = i1                  ! Index of VESSEL_WALL surface in the fluid code wall array
          surface(nsurface)%v(1,1) =  DBLE(wallpt(i2,20))  ! x coordinate of side vertex 1 (m)
          surface(nsurface)%v(2,1) =  DBLE(wallpt(i2,21))  ! y coordinate
          surface(nsurface)%v(3,1) = -1.0D+18              ! z (code currently assumes toroidal symmetry)
          surface(nsurface)%v(1,2) =  DBLE(wallpt(i1,20))  ! x coordinate of side vertex 2 (m)
          surface(nsurface)%v(2,2) =  DBLE(wallpt(i1,21))  ! y
          surface(nsurface)%v(3,2) =  1.0D+18              ! z
        ENDIF
      ENDDO

c...  Additional user specified wall surfaces (EIRASDAT array):
      DO i1 = 1, eirnasdat
        IF (eirasdat(i1,1) .EQ.-1.0.OR.
     .      eirasdat(i1,10).EQ. 0.0) CYCLE

        IF     (eirasdat(i1,1).EQ.7.0) THEN
c...      *** The 7.0 designation needs to change since not all the 7.0 info is being used ***
          nsurface = NewEireneSurface_06(VESSEL_WALL)

          surface(nsurface)%index(2) = NINT(eirasdat(i1,10))  ! Additional surface index

          surface(nsurface)%v(1,1) =  DBLE(eirasdat(i1,3))
          surface(nsurface)%v(2,1) =  DBLE(eirasdat(i1,4))
          surface(nsurface)%v(3,1) = -1.0D+18
          surface(nsurface)%v(1,2) =  DBLE(eirasdat(i1+1,3))
          surface(nsurface)%v(2,2) =  DBLE(eirasdat(i1+1,4))
          surface(nsurface)%v(3,2) =  1.0D+18

        ELSEIF (eirasdat(i1,1).EQ.11.0) THEN
c...      Tube (fun!):
          dangle = (eirasdat(i1,8) - eirasdat(i1,7)) / eirasdat(i1,9) 
          xcen = eirasdat(i1,2)
          ycen = eirasdat(i1,3)
          rad = eirasdat(i1,6)
          DO angle = eirasdat(i1,7), 0.9999*eirasdat(i1,8), dangle
            x1 = xcen + rad * COS(angle*PI/180.0)
            y1 = ycen + rad * SIN(angle*PI/180.0)
            x2 = xcen + rad * COS((angle+dangle)*PI/180.0)
            y2 = ycen + rad * SIN((angle+dangle)*PI/180.0)

            nsurface = NewEireneSurface_06(VESSEL_WALL)
            surface(nsurface)%index(2) = NINT(eirasdat(i1,10))  ! Additional surface index
            surface(nsurface)%v(1,1) = DBLE(x1)  ! *** change X1, etc. to REAL*8 ***
            surface(nsurface)%v(2,1) = DBLE(y1)
            surface(nsurface)%v(3,1) = DBLE(eirasdat(i1,4))
            surface(nsurface)%v(1,2) = DBLE(x2)
            surface(nsurface)%v(2,2) = DBLE(y2)
            surface(nsurface)%v(3,2) = DBLE(eirasdat(i1,5))
          ENDDO

        ELSEIF (eirasdat(i1,1).EQ.12.0) THEN
c...      Rectangular tube: 
          x1 = eirasdat(i1,2)
          y1 = eirasdat(i1,3)
          x2 = eirasdat(i1,4)
          y2 = eirasdat(i1,5)

          nsurface = NewEireneSurface_06(VESSEL_WALL)
          surface(nsurface)%index(2) = NINT(eirasdat(i1,10))  ! Additional surface index
          surface(nsurface)%v(1,1) = DBLE(x1)
          surface(nsurface)%v(2,1) = DBLE(y1)
          surface(nsurface)%v(3,1) = DBLE(eirasdat(i1,6))
          surface(nsurface)%v(1,2) = DBLE(x2)
          surface(nsurface)%v(2,2) = DBLE(y1)
          surface(nsurface)%v(3,2) = DBLE(eirasdat(i1,7))

          nsurface = NewEireneSurface_06(VESSEL_WALL)
          surface(nsurface)%index(2) = NINT(eirasdat(i1,10))  
          surface(nsurface)%v(1,1) = DBLE(x2)
          surface(nsurface)%v(2,1) = DBLE(y1)
          surface(nsurface)%v(3,1) = DBLE(eirasdat(i1,6))
          surface(nsurface)%v(1,2) = DBLE(x2)
          surface(nsurface)%v(2,2) = DBLE(y2)
          surface(nsurface)%v(3,2) = DBLE(eirasdat(i1,7))

          nsurface = NewEireneSurface_06(VESSEL_WALL)
          surface(nsurface)%index(2) = NINT(eirasdat(i1,10))  
          surface(nsurface)%v(1,1) = DBLE(x2)
          surface(nsurface)%v(2,1) = DBLE(y2)
          surface(nsurface)%v(3,1) = DBLE(eirasdat(i1,6))
          surface(nsurface)%v(1,2) = DBLE(x1)
          surface(nsurface)%v(2,2) = DBLE(y2)
          surface(nsurface)%v(3,2) = DBLE(eirasdat(i1,7))

          nsurface = NewEireneSurface_06(VESSEL_WALL)
          surface(nsurface)%index(2) = NINT(eirasdat(i1,10))  
          surface(nsurface)%v(1,1) = DBLE(x1)
          surface(nsurface)%v(2,1) = DBLE(y2)
          surface(nsurface)%v(3,1) = DBLE(eirasdat(i1,6))
          surface(nsurface)%v(1,2) = DBLE(x1)
          surface(nsurface)%v(2,2) = DBLE(y1)
          surface(nsurface)%v(3,2) = DBLE(eirasdat(i1,7))

        ELSEIF (eirasdat(i1,1).EQ.14.0) THEN
c...       
          DO i2 = i1, eirnasdat
            IF (i2.GT.i1.AND.eirasdat(i2,1).NE.0.0) EXIT
            IF (i2.EQ.i1) THEN
              x1 = eirasdat(i2,2)
              y1 = eirasdat(i2,3)
              z1 = eirasdat(i2,6)
              z2 = eirasdat(i2,7)
            ELSE
              x1 = x2
              y1 = y2
            ENDIF
            x2 = eirasdat(i2,4)
            y2 = eirasdat(i2,5)
            nsurface = NewEireneSurface_06(VESSEL_WALL)
            surface(nsurface)%index(2) = NINT(eirasdat(i1,10))  ! Additional surface index
            surface(nsurface)%v(1,1) =  DBLE(x1)
            surface(nsurface)%v(2,1) =  DBLE(y1)
            surface(nsurface)%v(3,1) =  DBLE(z1)
            surface(nsurface)%v(1,2) =  DBLE(x2)
            surface(nsurface)%v(2,2) =  DBLE(y2)
            surface(nsurface)%v(3,2) =  DBLE(z2)
          ENDDO

        ELSEIF (eirasdat(i1,1).EQ.15.0) THEN
c...      Holes in the triangle grid, so not really a surface, but no 
c         where else to put it:
          nsurface = NewEireneSurface_06(HOLE_IN_GRID)
          surface(nsurface)%index(2) = NINT(eirasdat(i1,10))  ! Additional surface index
          surface(nsurface)%index(3) = 1
          surface(nsurface)%v(1,1)   = DBLE(eirasdat(i1,2))
          surface(nsurface)%v(2,1)   = DBLE(eirasdat(i1,3))
        ELSE
        ENDIF
      ENDDO

c...  Additional user specified wall surfaces, new specification:
      DO i1 = 1, opt_eir%nadd
        assigned = .FALSE.
        SELECTCASE (opt_eir%add_type(i1))
c         --------------------------------------------------------------
          CASE (1,3)  ! Load line segment(s) and/or holes  from a file 
            fp = 99
            fname=TRIM(opt_eir%add_file    (i1))
            ftag =TRIM(opt_eir%add_file_tag(i1))
c            write(0,*) 'checking',opt_eir%add_type(i1),
c     .        TRIM(fname),TRIM(ftag)
            OPEN(UNIT=fp,FILE=TRIM(fname),ACCESS='SEQUENTIAL',
     .           STATUS='OLD',ERR=97)
            DO WHILE (.TRUE.)
              READ(fp,*) buffer
              IF (buffer(1:1).NE.'*') THEN
                READ(buffer,*) file_version                
c                WRITE(0,*) 'FILE VERSION!  ',file_version
                EXIT
              ENDIF
            ENDDO
            active = .FALSE.
            x2 = -999.0
            z1 = -1.0D+18
            z2 =  1.0D+18
            surface_type = VESSEL_WALL 
            IF (opt_eir%add_type(i1).EQ.3) THEN
              surface_type = HOLE_IN_GRID 
              hole_index   = 0
              file_version = 1.0
              x2 = 0.0
              y2 = 0.0
            ENDIF
            DO WHILE(osmGetLine(fp,buffer,ALL_LINES))
c              WRITE(0,*) 'BUFFER ',LEN_TRIM(buffer),'>'//
c     .                   TRIM(buffer)//'<'
              IF (buffer(1:1).EQ.'{'.AND.active) EXIT
              IF (active) THEN 
                SELECTCASE (NINT(file_version))
c                 --------------------------------------------------------
                  CASE (1)
                    READ(buffer,*,END=10) x1,y1
                    IF (x2.NE.-999.0) THEN 
c                      WRITE(0,*) 'processing...'
                      nsurface = NewEireneSurface_06(surface_type)
                      surface(nsurface)%index(2) = opt_eir%add_index(i1)
                      IF (surface_type.EQ.HOLE_IN_GRID) THEN
                        hole_index = hole_index + 1
                        surface(nsurface)%index(3) = hole_index
                      ENDIF
                      surface(nsurface)%v(1,1) = DBLE(x1)
                      surface(nsurface)%v(2,1) = DBLE(y1)
                      surface(nsurface)%v(3,1) = DBLE(z1)
                      surface(nsurface)%v(1,2) = DBLE(x2)
                      surface(nsurface)%v(2,2) = DBLE(y2)
                      surface(nsurface)%v(3,2) = DBLE(z2)
                    ENDIF
                    x2 = x1
                    y2 = y1
                    assigned(1) = .TRUE.
 10                 CONTINUE
c                 --------------------------------------------------------
                  CASE (2)
                    READ(buffer,*,END=20) x1,y1,x2,y2
                    nsurface = NewEireneSurface_06(surface_type)
c                    write(0,*) 'type',surface_type
c                    write(0,*) 'buff',nsurface,buffer(1:20)
                    surface(nsurface)%index(2) = opt_eir%add_index(i1)
                    surface(nsurface)%v(1,1) = DBLE(x1)
                    surface(nsurface)%v(2,1) = DBLE(y1)
                    surface(nsurface)%v(3,1) = DBLE(z1)
                    surface(nsurface)%v(1,2) = DBLE(x2)
                    surface(nsurface)%v(2,2) = DBLE(y2)
                    surface(nsurface)%v(3,2) = DBLE(z2)
                    assigned(1) = .TRUE.
 20                 CONTINUE
c                 ------------------------------------------------------
                  CASE DEFAULT
                    CALL ER('LoadEireneOption','Unknown additional '//
     .                      'line segment file version',*99)
c                 ------------------------------------------------------
                ENDSELECT
              ENDIF
              IF (buffer(1:1).EQ.'{'.AND..NOT.active.AND.
     .            osmCheckTag(buffer,ftag)) active = .TRUE.
            ENDDO
            CLOSE (fp)
c         --------------------------------------------------------------
          CASE (2)  ! A hole that was specified in the input file
            nsurface = NewEireneSurface_06(HOLE_IN_GRID)
            surface(nsurface)%index(2) = opt_eir%add_index(i1)
            surface(nsurface)%index(3) = 1
            surface(nsurface)%v(1,1)   = DBLE(opt_eir%add_holex(i1))
            surface(nsurface)%v(2,1)   = DBLE(opt_eir%add_holey(i1))
            assigned = .TRUE.
c         --------------------------------------------------------------
          CASE DEFAULT
            WRITE(0,*) 'add type= ',opt_eir%add_type(i1)
            CALL ER('LoadEireneOption','Unknown additional '//
     .              'surface type',*99)
        ENDSELECT
        IF (.NOT.assigned(1)) THEN
          write(0,*) '  ADD_TYPE = ',i1,opt_eir%add_type(i1)
          write(0,*) '  FNAME    = '//TRIM(fname)
          write(0,*) '  FTAG     = '//TRIM(ftag)
          CALL ER('DefineEireneSurfaces_06','Additional segment '//
     .            'wall tag not found',*99) 
        ENDIF
      ENDDO


      WRITE(eirfp,*) 'EIRENE SURFACE DEFINITION DATA - A:'
      DO i1 = 1, nsurface
        WRITE(eirfp,'(6I6,2X,A)')
     .    i1,
     .    surface(i1)%subtype,
     .    surface(i1)%index(1:2),
     .    surface(i1)%index(6),
     .    surface(i1)%iliin,                  
     .    TRIM(surface(i1)%surtxt)
      ENDDO


c     ------------------------------------------------------------------
c...  Over-ride default surface properties:
c     ------------------------------------------------------------------
      IF (opt_eir%sur_n.GT.0) THEN

        DO i1 = 1, nsurface
          DO i2 = 1, opt_eir%sur_n

            type = NINT(opt_eir%sur_type(i2))

c...        Identify which surface index to use for selecting surface
c           property specification:
            IF     (surface(i1)%type   .EQ.NON_DEFAULT_STANDARD  .AND.   ! Code not tested as I forgot that
     .              surface(i1)%subtype.EQ.MAGNETIC_GRID_BOUNDARY.AND.   ! I was using the template directly
     .              opt_eir%sur_type(i2).EQ.1.0) THEN
              index1 = surface(i1)%index(3)
            ELSEIF (surface(i1)%type   .EQ.NON_DEFAULT_STANDARD.AND.   
     .              surface(i1)%subtype.EQ.STRATUM             .AND.   
     .              opt_eir%sur_type(i2).EQ.1.1) THEN
              index1 = i1
            ELSEIF (surface(i1)%type.EQ.VESSEL_WALL.AND.
     .              opt_eir%sur_type(i2).EQ.2.0) THEN
              index1 = surface(i1)%index(1)
            ELSEIF (surface(i1)%type.EQ.VESSEL_WALL.AND.
     .              opt_eir%sur_type(i2).EQ.3.0) THEN
              index1 = surface(i1)%index(2)
            ELSE
              CYCLE
            ENDIF

            IF (index1.LE.0.OR.
     .          .NOT.CheckIndex(index1,0,opt_eir%sur_index(i2))) CYCLE
c            IF (.NOT.CheckIndex(index1,0,opt_eir%sur_index(i2))) CYCLE   ! changed 12/01/2016, SL
           
c            WRITE(0,*) 'THROUGH:',i1,i2,index1

            IF     (surface(i1)%type.EQ.NON_DEFAULT_STANDARD) THEN
              surface(i1)%surtxt = TRIM(surface(i1)%surtxt)//', '//
     .                             TRIM(opt_eir%sur_tag(i2))
            ELSEIF (surface(i1)%type.EQ.VESSEL_WALL         ) THEN
              surface(i1)%surtxt = TRIM(opt_eir%sur_tag(i2))
            ENDIF

            surface(i1)%iliin  = opt_eir%sur_iliin (i2)
            surface(i1)%ilside = opt_eir%sur_ilside(i2)
            surface(i1)%ilswch = opt_eir%sur_ilswch(i2)
            surface(i1)%recyct = opt_eir%sur_recyct(i2)

            surface(i1)%sector = opt_eir%sur_sector(i2)  ! Range of toroidal sectors for which the surface property is applied
            surface(i1)%hard   = opt_eir%sur_hard  (i2)  ! Force surface identity, i.e. no grouping (see below)

c...        Set surface sputtering parameters:
            IF (opt_eir%sur_ilspt(i2).NE.0) THEN
              tmp_ilspt = opt_eir%ilspt
              opt_eir%ilspt = 0
              IF (TRIM(opt_eir%sur_mat(i2)).EQ.'def') THEN
                opt_eir%ilspt = eirmat2
              ELSE
                IF (TRIM(opt_eir%sur_mat(i2)).EQ.'Mo') opt_eir%ilspt = 1
                IF (TRIM(opt_eir%sur_mat(i2)).EQ.'C' ) opt_eir%ilspt = 2
                IF (TRIM(opt_eir%sur_mat(i2)).EQ.'W' ) opt_eir%ilspt = 3
                IF (TRIM(opt_eir%sur_mat(i2)).EQ.'Be') opt_eir%ilspt = 4
                IF (TRIM(opt_eir%sur_mat(i2)).EQ.'Fe') opt_eir%ilspt = 5
                IF (opt_eir%ilspt.EQ.1) surface(i1)%material =  9642.0
                IF (opt_eir%ilspt.EQ.2) surface(i1)%material =  1206.0
                IF (opt_eir%ilspt.EQ.3) surface(i1)%material = 18474.0
                IF (opt_eir%ilspt.EQ.4) surface(i1)%material =   904.0
                IF (opt_eir%ilspt.EQ.5) surface(i1)%material =  5626.0
                surface(i1)%surtxt = TRIM(surface(i1)%surtxt)//', '//
     .                               TRIM(opt_eir%sur_mat(i2))
              ENDIF
              IF (opt_eir%ilspt.EQ.0) 
     .          CALL ER('DefineEireneSurfaces_06','Sputtered '//
     .                  'material not specified',*99)
              IF (tmp_ilspt.NE.0.AND.tmp_ilspt.NE.opt_eir%ilspt)    ! Can only have one type of material being sputtered
     .          CALL ER('DefineEireneSurfaces_06','Sputtered '//    ! in EIRENE at the moment
     .                  'material over-specified (not allowed)',*99)

              surface(i1)%ilspt  = opt_eir%sur_ilspt(i2)
              surface(i1)%isrs   = 2   ! Species index of sputtered atom
              surface(i1)%recycs = 1.0
              surface(i1)%recycc = 1.0
              IF (first_message) THEN
                WRITE(0,*) '*** SPUTTERING ON IN EIRENE ***',
     .                     surface(i1)%ilspt
                first_message = .FALSE.
              ENDIF
            ENDIF

c...        Set local temperature:
            IF (opt_eir%sur_temp(i2).GT.0) 
     .        surface(i1)%ewall = -opt_eir%sur_temp(i2) * 1.38E-23 / ECH

            surface(i1)%transp1 = opt_eir%sur_tr1(i2)  ! Surface transparency (default = 0.0, i.e. fully reflecting)
            surface(i1)%transp2 = opt_eir%sur_tr2(i2)

c...        Setup forced remapping of conformal surfaces, i.e. when the poloidal
c           boundary of a ring is coincident with the vessel wall specification, as
c           can happen with extended grids:
            IF (opt_eir%sur_remap(i2).GT.0) THEN
              IF (surface(i1)%type   .EQ.NON_DEFAULT_STANDARD  .AND.   
     .            surface(i1)%subtype.EQ.MAGNETIC_GRID_BOUNDARY) THEN
c               ***MASSIVE HACK : CONFORMAL SURFACE REMAPPING, SEE PROCESSTRIANGLES_06
                surface(i1)%index(10) = opt_eir%sur_remap(i2)
                WRITE(0,*) '*** SUPER HACK ACTIVATED! ***',i1
              ENDIF
            ENDIF

          ENDDO
        ENDDO

      ELSE

c...    Old method:
        DO i1 = 1, nsurface
        
          DO i2 = 1, eirnspdat
            type   = NINT(eirspdat(i2,1))
            index1 = NINT(eirspdat(i2,2))
            IF (eirspdat(MIN(i2+1,eirnspdat),1).EQ.0.0) THEN
              index2 = NINT(eirspdat(i2+1,2))
            ELSE
              index2 = index1
            ENDIF
            IF (type.EQ.0) CYCLE
        
            IF ((surface(i1)%type   .EQ.NON_DEFAULT_STANDARD  .AND.   ! Code not tested as I forgot that
     .           surface(i1)%subtype.EQ.MAGNETIC_GRID_BOUNDARY.AND.   ! I was using the template directly
     .             eirspdat(i2,1).EQ.1.0.AND.
     .             surface(i1)%index(3).GE.index1.AND.
     .             surface(i1)%index(3).LE.index2).OR.  ! *whew!*
     .          (surface(i1)%type   .EQ.NON_DEFAULT_STANDARD.AND.   
     .           surface(i1)%subtype.EQ.STRATUM             .AND.   
     .             eirspdat(i2,1).EQ.1.1.AND.
     .             i1.GE.index1.AND.
     .             i1.LE.index2).OR.  ! *whew!*
     .          (surface(i1)%type.EQ.VESSEL_WALL.AND.
     .           ((type.EQ.2.AND.
     .             surface(i1)%index(1).GE.index1.AND.
     .             surface(i1)%index(1).LE.index2).OR.
     .            (type.EQ.3.AND.
     .             surface(i1)%index(2).GE.index1.AND.
     .             surface(i1)%index(2).LE.index2)))) THEN
        
              surface(i1)%iliin  = NINT(eirspdat(i2,3))
              surface(i1)%ilside = NINT(eirspdat(i2,4))
              surface(i1)%ilswch = NINT(eirspdat(i2,5))
              surface(i1)%recyct = eirspdat(i2,8)

c             *** HACK FOR SPUTTERING AND CONFORMAL SURFACE REMAPPING***
              IF (eirspdat(i2,10).LT.0.0) THEN   ! (,10) is remapped from (,9) in unstructured_input_com.f
                IF (surface(i1)%type   .EQ.NON_DEFAULT_STANDARD  .AND.   
     .              surface(i1)%subtype.EQ.MAGNETIC_GRID_BOUNDARY) THEN
c                 *** MASSIVE HACK : SURFACE REMAPPING, SEE PROCESSTRIANGLES_06
                  surface(i1)%index(10) = -NINT(eirspdat(i2,10))
                  WRITE(0,*) '*** SUPER HACK ACTIVATED! ***',i1
c                surface(i1)%recyct = 1.0
                ELSE
                  opt_eir%ilspt = -NINT(eirspdat(i2,10)) 
                  IF (opt_eir%ilspt.EQ.1) opt_eir%ilspt = 2  ! This remapping is necessary because of the 
                  IF (opt_eir%ilspt.EQ.2) opt_eir%ilspt = 5  ! original sputtering hack I forced into the
                  IF (opt_eir%ilspt.EQ.3) opt_eir%ilspt = 3  ! old method of specifying the surface 
                  IF (opt_eir%ilspt.EQ.4) opt_eir%ilspt = 4  ! properties.
                  surface(i1)%ilspt  = 2  
                  surface(i1)%isrs   = 2   ! Species index of sputtered atom
                  surface(i1)%recycs = 1.0
                  surface(i1)%recycc = 1.0
                  WRITE(0,*) '*** SPUTTERING ON IN EIRENE ***',
     .                       surface(i1)%ilspt
                ENDIF
              ENDIF   
            ENDIF
            ENDDO
        ENDDO
      ENDIF
c
c     If the above surfaces are defined correctly, then the following
c     should be the same for all plasma codes.
c
c...  Count the number of non-default standard surfaces already 
c     defined (magnetic grid boundaries specified above):
      sur1 = 0
      DO i1 = 1, nsurface
        IF (surface(i1)%type.EQ.NON_DEFAULT_STANDARD) sur1 = sur1 + 1
      ENDDO
      sur1 = sur1 + 1
      sur2 = sur1 - 1
c...  Build non-default standard surfaces to represent vessel wall surfaces when
c     writing block 3 of the the EIRENE06 input file.  A different surface is
c     defined each time a different surface reflection model and/or surface 
c     property (material, temperature, etc.) is identified.  At present, only the
c     TYPE.EQ.NON_DEFAULT_STANDARD surfaces in the list of surfaces are used
c     when writing block 3a in the EIRENE input file (and block 3b currently has
c     zero surfaces):
      ntmp = nsurface
c...  Set the default wall surface:
      nsurface = NewEireneSurface_06(NON_DEFAULT_STANDARD)
      default_surface = nsurface
      sur2 = sur2 + 1
      surface(nsurface)%subtype  = ADDITIONAL
      surface(nsurface)%index(5) = sur2
      surface(nsurface)%reflect  = LOCAL
      surface(nsurface)%ewall    = -wtemp * 1.38E-23 / ECH
      surface(nsurface)%material = wmater
      surface(nsurface)%surtxt   = '* default wall surface (DIVIMP)'
c...  Scan through vessel wall surfaces and assign them to a non-default standard
c     EIRENE surface, creating new non-default surfaces as required:
      DO i1 = 1, ntmp
        IF (surface(i1)%type.NE.VESSEL_WALL) CYCLE
c...    Surface matching criteria:
        iliin   = surface(i1)%iliin
        ilside  = surface(i1)%ilside
        ilswch  = surface(i1)%ilswch
        ilspt   = surface(i1)%ilspt
        mater   = surface(i1)%material
        ewall   = surface(i1)%ewall
        transp1 = surface(i1)%transp1
        transp2 = surface(i1)%transp2
        recycf  = surface(i1)%recycf
        recyct  = surface(i1)%recyct
        isrs    = surface(i1)%isrs 
        recycs  = surface(i1)%recycs
        recycc  = surface(i1)%recycc
c...    Check if a non-default standard surface has already been defined
c       that matches the above citeria:
        sur3 = 0

c        WRITE(0,*) 'sur:',i1
        DO i2 = ntmp+1, nsurface
 
c          WRITE(0,'(I5,11L2)')  i2,
c     .     surface(i2)%iliin   .EQ.iliin, 
c     .     surface(i2)%ilside  .EQ.ilside,
c     .     surface(i2)%ilswch  .EQ.ilswch,
c     .     surface(i2)%ilspt   .EQ.ilspt, 
c     .     surface(i2)%material.EQ.mater, 
c     .     surface(i2)%ewall   .EQ.ewall, 
c     .     surface(i2)%recycf  .EQ.recycf,
c     .     surface(i2)%recyct  .EQ.recyct,
c     .     surface(i2)%isrs    .EQ.isrs,  
c     .     surface(i2)%recycs  .EQ.recycs,
c     .     surface(i2)%recycc  .EQ.recycc

c         WRITE(0,*) 
c     .     surface(i2)%ewall,ewall 

          IF (surface(i2)%iliin   .EQ.iliin  .AND.
     .        surface(i2)%ilside  .EQ.ilside .AND.
     .        surface(i2)%ilswch  .EQ.ilswch .AND.
     .        surface(i2)%ilspt   .EQ.ilspt  .AND.
     .        surface(i2)%material.EQ.mater  .AND.
     .        surface(i2)%ewall   .EQ.ewall  .AND.
     .        surface(i2)%transp1 .EQ.transp1.AND.
     .        surface(i2)%transp2 .EQ.transp2.AND.
     .        surface(i2)%recycf  .EQ.recycf .AND.
     .        surface(i2)%recyct  .EQ.recyct .AND.
     .        surface(i2)%isrs    .EQ.isrs   .AND.
     .        surface(i2)%recycs  .EQ.recycs .AND.
     .        surface(i2)%recycc  .EQ.recycc .AND.
     .        surface(i2)%hard    .EQ.0) sur3 = i2
        ENDDO                    
        IF (sur3.EQ.0.OR.surface(i1)%hard.NE.0) THEN
c...      Existing surface matching the criteria for the current surface 
c         not found, or separate surface represenation forced, so add a surface:
          sur2 = sur2 + 1
          nsurface = NewEireneSurface_06(NON_DEFAULT_STANDARD)
          surface(nsurface)%subtype  = ADDITIONAL
          surface(nsurface)%iliin    = surface(i1)%iliin
          surface(nsurface)%ilside   = surface(i1)%ilside
          surface(nsurface)%ilswch   = surface(i1)%ilswch
          surface(nsurface)%ilspt    = surface(i1)%ilspt 
          surface(nsurface)%index(5) = sur2
          surface(nsurface)%reflect  = LOCAL
          surface(nsurface)%material = surface(i1)%material
          surface(nsurface)%ewall    = surface(i1)%ewall
          surface(nsurface)%transp1  = surface(i1)%transp1
          surface(nsurface)%transp2  = surface(i1)%transp2
          surface(nsurface)%recyct   = surface(i1)%recyct
          surface(nsurface)%recycf   = surface(i1)%recycf
          surface(nsurface)%isrs     = surface(i1)%isrs  
          surface(nsurface)%recycs   = surface(i1)%recycs
          surface(nsurface)%recycc   = surface(i1)%recycc
          surface(nsurface)%hard     = surface(i1)%hard  
          surface(nsurface)%sector   = surface(i1)%sector
          SELECTCASE (iliin)
            CASE (-2)
              surface(nsurface)%surtxt = 
     .          '* transparent switching surface (DIVIMP)'
            CASE (0)
              surface(nsurface)%surtxt = 
     .          '* transparent non-switching surface (DIVIMP)'
            CASE (1)
              surface(nsurface)%surtxt = '* wall surface (DIVIMP)'      
            CASE (2)                
              surface(nsurface)%surtxt = '* pumping surface (DIVIMP)'      
            CASE (3)                
              surface(nsurface)%surtxt ='* specular reflection (DIVIMP)'      
            CASE DEFAULT
              CALL ER('DefineEireneSurfaces','Invalid ILIIN',*99)
          ENDSELECT
          IF (surface(nsurface)%hard.NE.0) 
     .      surface(nsurface)%surtxt = 
     .        TRIM(surface(nsurface)%surtxt)//', '//
     .        TRIM(surface(i1      )%surtxt)
c...      Assign the wall surface to the non-default standard surface:
          surface(i1)%index(3) = sur2
        ELSE
c...      Assign the current wall surface to the identified non-default
c         standard surface that already exists:
          surface(i1)%index(3) = surface(sur3)%index(5) 
        ENDIF
      ENDDO

c...  Add a surface for toroidal boundaries of tetrahedral grids -- should eventually make
c     periodic in EIRENE:
      IF (tetrahedrons) THEN
        nsurface = NewEireneSurface_06(NON_DEFAULT_STANDARD)
        surface(nsurface)%subtype  = ADDITIONAL
        surface(nsurface)%index(1) = -1
        surface(nsurface)%reflect  = LOCAL
        surface(nsurface)%ewall    = -wtemp * 1.38E-23 / ECH
        surface(nsurface)%material = wmater
        surface(nsurface)%iliin    = opt_eir%tet_iliin
        surface(nsurface)%ilcol    = 5         
        surface(nsurface)%surtxt   = '* tetrahedron dump surface (DIV)'
      ELSE
        nsurface = NewEireneSurface_06(NON_DEFAULT_STANDARD)
        surface(nsurface)%subtype  = ADDITIONAL
        surface(nsurface)%index(1) = -1
        surface(nsurface)%reflect  = LOCAL
        surface(nsurface)%ewall    = -wtemp * 1.38E-23 / ECH
        surface(nsurface)%material = wmater
        surface(nsurface)%iliin    = 1
        surface(nsurface)%ilcol    = 5         
        surface(nsurface)%surtxt   = '* triangle dump surface (DIVIMP)'
      ENDIF

c...  Assign block 3a surface index to non-default standard surfaces:
      i2 = 0
      DO i1 = 1, nsurface
        IF (surface(i1)%type.EQ.NON_DEFAULT_STANDARD) THEN
          i2 = i2 + 1
          surface(i1)%num = i2
        ENDIF
      ENDDO

c...  Output:
      WRITE(eirfp,*) 'EIRENE SURFACE DEFINITION DATA - FINAL:'
      DO i1 = 1, nsurface
        WRITE(eirfp,'(6I6,2X,A)')
     .    i1,
     .    surface(i1)%subtype,
     .    surface(i1)%index(1:2),
     .    surface(i1)%index(6),
     .    surface(i1)%iliin,                  
     .    TRIM(surface(i1)%surtxt)
      ENDDO

      RETURN
 97   CALL ER('DefineEireneSurfaces_06','Unable to open file',*98)
 98   WRITE(0,*) '  FILE= '//TRIM(fname)
      STOP
 99   WRITE(0,*) '  ILIIN: ',iliin
      STOP
      END
c
c ======================================================================
c
c subroutine: ProcessFluidGrid
c
c   Cell organization in OSM:
c
c        TARGET 2 (=IKHI)
c
c        *-------------*
c        |             |
c        | IK=NKS(IR)  |
c        |             |
c 
c              ...
c 
c        |             |
c        |    IK=2     |
c        |             |
c  VER 3 *-------------* VERTEX 4
c        |   SIDE 3    |
c        |S           S|
c        |I           I|
c        |D           D|
c        |E   IK=1    E|
c        |             |
c        |2           4|
c        |   SIDE 1    |
c        *-------------*
c   VERTEX 2        VERTEX 1
c
c        TARGET 1 (=IKLO)
c
c
      SUBROUTINE ProcessFluidGrid_06
      USE mod_eirene06_parameters
      USE mod_eirene06_locals
      USE mod_eirene06
      USE mod_grid_divimp
      use mod_params
      use mod_comtor
      use mod_cgeom
      use mod_slcom
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'comtor'
c     INCLUDE 'cgeom'
c     INCLUDE 'slcom'

      REAL       TOL
      PARAMETER (TOL=1.0E-07)

      REAL GetMach,GetJsat,GetFlux 

      INTEGER ike,ik,ir,id,i1,i2,i3,in,it,region,ik1,ir1,path,fp,code,
     .        iw,is1,is2
      REAL    fact,Bfrac,e_pot(0:MAXNKS+1,MAXNRS),frac,maxepot,area
      REAL*8  Bx,By,Bz,beta,brat,deltax,deltay,x(3),y(3),
     .        tri_r1,tri_z1,tri_r2,tri_z2

      INTEGER :: be_n(0:10),be_ir_n,be_ir(10)
c      INTEGER :: wall_n
      REAL*8, PARAMETER :: DTOL = 1.0D-06

      TYPE :: bnd_element
        INTEGER :: ik,ir,c,sideindex(5)
        REAL*8  :: r(2),z(2)
      ENDTYPE bnd_element
      TYPE(bnd_element), ALLOCATABLE :: be(:)

      CHARACTER command*512

      INTEGER, PARAMETER :: MAX_BE_N = 10000
c      INTEGER, PARAMETER :: MAX_WALL_N = 10000

c      WRITE(0,*) 'PROCESSING MAGNETIC GRID'

      IF (fluid_grid.EQ.2.OR.fluid_grid.EQ.3) THEN  ! fluid triangles
c       -----------------------------------------------------------------
c...    Triangular fluid grid proxy:

c...    Trace the outer boundary of the fluid grid:
        be_n = 0
        ALLOCATE(be(0:MAX_BE_N))

c        ALLOCATE(wall_r (2,0:MAX_WALL_N))
c        ALLOCATE(wall_z (2,0:MAX_WALL_N))
c        ALLOCATE(wall_ik(  0:MAX_WALL_N))
c        ALLOCATE(wall_ir(  0:MAX_WALL_N))
c        ALLOCATE(wall_c (  0:MAX_WALL_N))

         be(:)%ik = 0
         be(:)%ir = 0
         be(:)%c  = 0
         DO i1 = 1, MAX_BE_N
           be(i1)%sideindex(1:5) = 0
         ENDDO
c        wall_ik = 0
c        wall_ir = 0

c...    Build list of segments:
        DO ir = irsep, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE

          IF (be_n(0).GE.MAX_BE_N-2) 
     .      CALL ER('ProcessFluidGrid_06','MAX_BE_N exceeded (1)',*99) 

          be_n = be_n(0) + 1
          be(be_n(0))%ik = 1
          be(be_n(0))%ir = ir
          be(be_n(0))%sideindex(2) = IKLO
          be(be_n(0))%sideindex(5) = nimindex(idds(ir,2))
          id = korpg(1,ir)
          IF (ALLOCATED(d_rvertp)) THEN
            be(be_n(0))%r(1) = d_rvertp(1,id)       ! x (m)
            be(be_n(0))%z(1) = d_zvertp(1,id)       ! y
            be(be_n(0))%r(2) = d_rvertp(2,id)       ! x (m)
            be(be_n(0))%z(2) = d_zvertp(2,id)       ! y
            be_n(0) = be_n(0) + 1
            id = korpg(nks(ir),ir)
            be(be_n(0))%r(1) = d_rvertp(3,id)
            be(be_n(0))%z(1) = d_zvertp(3,id)
            be(be_n(0))%r(2) = d_rvertp(4,id)
            be(be_n(0))%z(2) = d_zvertp(4,id)
          ELSE
            be(be_n(0))%r(1) = DBLE(rvertp(1,id))           
            be(be_n(0))%z(1) = DBLE(zvertp(1,id))
            be(be_n(0))%r(2) = DBLE(rvertp(2,id))           
            be(be_n(0))%z(2) = DBLE(zvertp(2,id))
            be_n(0) = be_n(0) + 1
            id = korpg(nks(ir),ir)
            be(be_n(0))%r(1) = DBLE(rvertp(3,id))           
            be(be_n(0))%z(1) = DBLE(zvertp(3,id))
            be(be_n(0))%r(2) = DBLE(rvertp(4,id))           
            be(be_n(0))%z(2) = DBLE(zvertp(4,id))
          ENDIF
          be(be_n(0))%ik = nks(ir)
          be(be_n(0))%ir = ir
          be(be_n(0))%sideindex(2) = IKHI
          be(be_n(0))%sideindex(5) = nimindex(idds(ir,1))
        ENDDO 

        DO in = 1, 2
          IF (in.EQ.1) ir1 = IRWALL
          IF (in.EQ.2) ir1 = IRTRAP
          IF (ir1.GT.nrs) EXIT

          IF (be_n(0).GE.MAX_BE_N-nks(ir1))
     .      CALL ER('ProcessFluidGrid_06','MAX_BE_N exceeded (2)',*99) 
 
          DO ik1 = 1, nks(ir1)
            IF (in.EQ.1) THEN
              ik = ikins(ik1,ir1)
              ir = irins(ik1,ir1)
            ELSE
              ik = ikouts(ik1,ir1)
              ir = irouts(ik1,ir1)
            ENDIF
c            write(0,*) ik,ir,irouts(ik,ir).EQ.ir1
            id = korpg(ik,ir)
            be_n(0) = be_n(0) + 1
            be(be_n(0))%ik = ik
            be(be_n(0))%ir = ir
            IF (irouts(ik,ir).EQ.ir1) THEN
              be(be_n(0))%sideindex(1) = 23
              IF (ALLOCATED(d_rvertp)) THEN
                be(be_n(0))%r(1) = d_rvertp(2,id)
                be(be_n(0))%z(1) = d_zvertp(2,id)
                be(be_n(0))%r(2) = d_rvertp(3,id)
                be(be_n(0))%z(2) = d_zvertp(3,id)
              ELSE
                be(be_n(0))%r(1) = DBLE(rvertp(2,id))           
                be(be_n(0))%z(1) = DBLE(zvertp(2,id))
                be(be_n(0))%r(2) = DBLE(rvertp(3,id))           
                be(be_n(0))%z(2) = DBLE(zvertp(3,id))
              ENDIF     
            ELSE
              be(be_n(0))%sideindex(1) = 14
              IF (ALLOCATED(d_rvertp)) THEN
                be(be_n(0))%r(1) = d_rvertp(4,id)
                be(be_n(0))%z(1) = d_zvertp(4,id)
                be(be_n(0))%r(2) = d_rvertp(1,id)
                be(be_n(0))%z(2) = d_zvertp(1,id)
              ELSE
                be(be_n(0))%r(1) = DBLE(rvertp(4,id))           
                be(be_n(0))%z(1) = DBLE(zvertp(4,id))
                be(be_n(0))%r(2) = DBLE(rvertp(1,id))           
                be(be_n(0))%z(2) = DBLE(zvertp(1,id))
              ENDIF     
            ENDIF
          ENDDO
        ENDDO

c....   Sort the segments:
        WRITE(88,'(A,2(2F14.10,2X))') 
     .    'sorting',be(1)%r(1),be(1)%z(1),
     .              be(1)%r(2),be(1)%z(2)

        DO i1 = 1, be_n(0)-1
          DO i2 = i1+1, be_n(0)
            IF (DABS(be(i1)%r(2)-be(i2)%r(1)).LT.DTOL.AND.
     .          DABS(be(i1)%z(2)-be(i2)%z(1)).LT.DTOL) THEN
              IF (i2.NE.i1+1) THEN
c               Swap:
                be(0   ) = be(i1+1)
                be(i1+1) = be(i2  )
                be(i2  ) = be(0   )
c                wall_r (:,0   ) = wall_r(:,i1+1)
c                wall_z (:,0   ) = wall_z(:,i1+1)
c                wall_ik(0   ) = wall_ik(i1+1)
c                wall_ir(0   ) = wall_ir(i1+1)
c                wall_r(:,i1+1) = wall_r(:,i2  )
c                wall_z(:,i1+1) = wall_z(:,i2  )
c                wall_r(:,i2  ) = wall_r(:,0   )
c                wall_z(:,i2  ) = wall_z(:,0   )                    
              ENDIF
              WRITE(88,'(A,2(2F14.10,2X))') 
     .          'sorting',be(i1+1)%r(1),be(i1+1)%z(1),
     .                    be(i1+1)%r(2),be(i1+1)%z(2)
              EXIT
            ENDIF
          ENDDO
          IF (i2.EQ.be_n(0)+1) 
     .      CALL ER('ProcessFluidGrid_06','Wall gap detected',*99) 
        ENDDO

        IF (cgridopt.EQ.LINEAR_GRID) THEN
          be_n(0) = be_n(0) + 1
          be(be_n(0))%r(1) = be(be_n(0)-1)%r(2)
          be(be_n(0))%z(1) = be(be_n(0)-1)%z(2)
          be(be_n(0))%r(2) = be(1     )%r(1)
          be(be_n(0))%z(2) = be(1     )%z(1)
        ENDIF

        DO i1 = 1, be_n(0)
          WRITE(88,'(A,2(2F14.10,2X))') 
     .      'sort list',be(i1)%r(1),be(i1)%z(1),
     .                  be(i1)%r(2),be(i1)%z(2)
        ENDDO

        be_n(1) = be_n(0)

        IF (fluid_grid.EQ.3) THEN
c...      Add the separatrix:      

          be_ir_n = 5
          be_ir(1:be_ir_n) = (/2,irsep-9,irsep,irsep+14,47/)

          DO i1 = 1, be_ir_n
            ir = be_ir(i1)
            ike = nks(ir)-1
            IF (ir.LT.irsep) ike = ike - 1
            DO ik = 1, ike
              be_n(0) = be_n(0) + 1
              be(be_n(0))%r(1) = DBLE(rvertp(4,korpg(ik,ir)))
              be(be_n(0))%z(1) = DBLE(zvertp(4,korpg(ik,ir)))
c              be(be_n(0))%r(1) = DBLE(rs(ik,ir))  
c              be(be_n(0))%z(1) = DBLE(zs(ik,ir))
            ENDDO
            be_n(i1+1) = be_n(0)          
          ENDDO

c          DO ik = 1, nks(irsep)-1
c            be_n(0) = be_n(0) + 1
c            be(be_n(0))%r(1) = DBLE(rvertp(4,korpg(ik,irsep)))
c            be(be_n(0))%z(1) = DBLE(zvertp(4,korpg(ik,irsep)))
cc            be(be_n(0))%r(1) = DBLE(rs(ik,irsep))  
cc            be(be_n(0))%z(1) = DBLE(zs(ik,irsep))
c          ENDDO
c          be_n(2) = be_n(0)

c          DO ik = 1, nks(irsep-15)-2
c            be_n(0) = be_n(0) + 1
c            be(be_n(0))%r(1) = DBLE(rvertp(4,korpg(ik,irsep-15)))
c            be(be_n(0))%z(1) = DBLE(zvertp(4,korpg(ik,irsep-15)))
c          ENDDO
c          be_n(3) = be_n(0)

c          DO ik = 1, nks(irsep-9)-2
c            be_n(0) = be_n(0) + 1
c            be(be_n(0))%r(1) = DBLE(rvertp(4,korpg(ik,irsep-9)))
c            be(be_n(0))%z(1) = DBLE(zvertp(4,korpg(ik,irsep-9)))
c          ENDDO
c          be_n(4) = be_n(0)

c          DO ik = 1, nks(irsep+10)-1
c            be_n(0) = be_n(0) + 1
c            be(be_n(0))%r(1) = DBLE(rvertp(4,korpg(ik,irsep+10)))
c            be(be_n(0))%z(1) = DBLE(zvertp(4,korpg(ik,irsep+10)))
c          ENDDO
c          be_n(5) = be_n(0)

c          DO ik = 1, nks(irsep+20)-1
c            be_n(0) = be_n(0) + 1
c            be(be_n(0))%r(1) = DBLE(rvertp(4,korpg(ik,irsep+20)))
c            be(be_n(0))%z(1) = DBLE(zvertp(4,korpg(ik,irsep+20)))
c          ENDDO
c          be_n(6) = be_n(0)

        ENDIF


c...    Call TRIANGLE:
        fp = 99      
        OPEN(UNIT=fp,FILE='triangle.poly',ACCESS='SEQUENTIAL',
     .       STATUS='REPLACE',ERR=99)      
        WRITE(fp,*) be_n(0),2,1,0
        DO i1 = 1, be_n(0)
          WRITE(fp,'(I6,2F19.14)') i1,be(i1)%r(1),be(i1)%z(1)
        ENDDO

        IF (fluid_grid.EQ.3) THEN

          WRITE(fp,*) be_n(be_ir_n+1)-be_ir_n,0  ! *** here, two times ***
c          WRITE(fp,*) be_n(6)-1-1-1-1-1,0  ! *** here, two times ***

          DO i2 = 1, be_n(1)
            IF (i2.EQ.be_n(1)) THEN
              WRITE(fp,'(4I6)') i2,be_n(1),1   ,0
            ELSE
              WRITE(fp,'(4I6)') i2,i2     ,i2+1,1
            ENDIF
          ENDDO
          DO i3 = 1, be_ir_n   ! *** here ***
            DO i2 = be_n(i3)+1, be_n(i3+1)-1
              IF (i2.EQ.be_n(i3+1)-1) THEN
                WRITE(fp,'(4I6)') i2,i2   ,i2+1   ,0
              ELSE
                WRITE(fp,'(4I6)') i2,i2   ,i2+1   ,1
              ENDIF
            ENDDO
          ENDDO

        ELSE

          WRITE(fp,*) be_n(1),0  ! *** here, two times ***

          DO i2 = 1, be_n(1)
            IF (i2.EQ.be_n(1)) THEN
              WRITE(fp,'(4I6)') i2,be_n(1),1   ,0
            ELSE
              WRITE(fp,'(4I6)') i2,i2     ,i2+1,1
            ENDIF
          ENDDO

        ENDIF

        WRITE(fp,*) 0 ! nhole
c        DO i2 = 1, nhole
c          WRITE(fp,'(I6,2F12.7)') i2,xhole(i2),yhole(i2)
c        ENDDO
        CLOSE (fp)
       
c.....  Call triangle:
        area = 1.0
        WRITE(command,10) 'triangle -p -q -a',area,
     .                    ' -Y triangle.poly>tmp'
 10     FORMAT(A,F10.8,A)
        WRITE(eirfp,*) 'COMMAND: >'//command(1:LEN_TRIM(command))//'<'
        WRITE(0    ,*) 'COMMAND: >'//command(1:LEN_TRIM(command))//'<'
        CALL CIssue(command(1:LEN_TRIM(command)),code)
        WRITE(eirfp,*) 'RETURN_CODE:',code
      
        write(0,*) 'ntri=',ntri

        CALL ReadPolyFile_06(0,1.0E+12,10.0,10.0)

        write(0,*) 'ntri=',ntri

c...    Map the triangles to fluid cells: 
        DO it = 1, ntri
          DO is1 = 1, 3
            is2 = is1 + 1
            IF (is2.EQ.4) is2 = 1

            tri_r1 = ver(tri(it)%ver(is1),1)
            tri_z1 = ver(tri(it)%ver(is1),2)
            tri_r2 = ver(tri(it)%ver(is2),1)
            tri_z2 = ver(tri(it)%ver(is2),2)

            DO iw = 1, be_n(0)          

              IF ((DABS(be(iw)%r(1)-tri_r1).LT.DTOL.AND.
     .             DABS(be(iw)%z(1)-tri_z1).LT.DTOL.AND.
     .             DABS(be(iw)%r(2)-tri_r2).LT.DTOL.AND.
     .             DABS(be(iw)%z(2)-tri_z2).LT.DTOL).OR.
     .            (DABS(be(iw)%r(2)-tri_r1).LT.DTOL.AND.
     .             DABS(be(iw)%z(2)-tri_z1).LT.DTOL.AND.
     .             DABS(be(iw)%r(1)-tri_r2).LT.DTOL.AND.
     .             DABS(be(iw)%z(1)-tri_z2).LT.DTOL)) THEN

                be(iw)%c = be(iw)%c + 1

                tri(it)%index(IND_IK) = be(iw)%ik
                tri(it)%index(IND_IR) = be(iw)%ir

                tri(it)%sideindex(1,is1) = be(iw)%sideindex(1)
                tri(it)%sideindex(2,is1) = be(iw)%sideindex(2)
                tri(it)%sideindex(5,is1) = be(iw)%sideindex(5)

              ENDIF
            ENDDO

          ENDDO
        ENDDO

        DO iw = 1, be_n(0)          
          IF (be(iw)%ik.NE.0.AND.be(iw)%c.NE.1) THEN
            WRITE(0,*) 'SHIT',iw,be(iw)%c
          ENDIF
        ENDDO

c...    Load cell geometry and volume quantities for the magnetic/fluid grid:
c        fact = qtim * qtim * emi / crmi
        fact = 1.0
        ncell = 0
        DO it = 1, ntri
          ik = tri(it)%index(IND_IK)
          ir = tri(it)%index(IND_IR)

          ncell = ncell + 1
          cell(ncell)%index     = 0
          cell(ncell)%sideindex = 0
          cell(ncell)%type = 1
c...      Cell indices on the magnetic grid: 
          cell(ncell)%index(1) = ik                                ! Cell index
          cell(ncell)%index(2) = ir                                ! Ring index
c...      Radial cell surfaces:
          cell(ncell)%sideindex(1,1:3) = tri(it)%sideindex(1,1:3)  ! Cell side index
          cell(ncell)%sideindex(2,1:3) = tri(it)%sideindex(2,1:3)  ! Target index 
          cell(ncell)%sideindex(5,1:3) = tri(it)%sideindex(5,1:3)  ! Map to DIVIMP wall array
c...      Cell vertices (3 vertices assumed at present):            
          cell(ncell)%r(1:3) = ver(tri(it)%ver(1:3),1)
          cell(ncell)%z(1:3) = ver(tri(it)%ver(1:3),2)
          IF (ik.NE.0) THEN
c...        B-field components (approximate):         
            Bz = 0.0
            By = 0.0
            Bx = 0.0
c...        CBPHI is the on-axis B-field value specified in the OSM input file:
            Bfrac = cbphi * r0 / rs(ik,ir)                         ! Rough scaling of B-field
c...        Plasma quantities:
            cell(ncell)%plasma(1) = ktebs(ik,ir)                   ! Te (eV)
            cell(ncell)%plasma(2) = ktibs(ik,ir)                   ! Ti (eV)
            cell(ncell)%plasma(3) = knbs (ik,ir)                   ! ni (eV) (ne=ni assumed at present)
            cell(ncell)%plasma(4) = SNGL(Bx) * kvhs(ik,ir)         ! vx (m-1 s-1)
            cell(ncell)%plasma(5) = SNGL(By) * kvhs(ik,ir)         ! vy
            cell(ncell)%plasma(6) = SNGL(Bz) * kvhs(ik,ir)         ! vz
c...        E&M quantities:
            cell(ncell)%bfield(1) = SNGL(Bx) * Bfrac               ! Bx (Tesla) (normalized on Eirene side)
            cell(ncell)%bfield(2) = SNGL(By) * Bfrac               ! By 
            cell(ncell)%bfield(3) = SNGL(Bz) * Bfrac               ! Bz 
            cell(ncell)%bfield(4) = brat                           ! Bratio 
            cell(ncell)%efield(1) = SNGL(Bx) * kes(ik,ir) / fact   ! Ex (not required by EIRENE)
            cell(ncell)%efield(2) = SNGL(By) * kes(ik,ir) / fact   ! Ey
            cell(ncell)%efield(3) = SNGL(Bz) * kes(ik,ir) / fact   ! Ez
            
            cell(ncell)%e_pot = e_pot(ik,ir)                       ! Electric potential (estimate)
          ELSE
            cell(ncell)%plasma = 0.0
            cell(ncell)%bfield = 0.0
            cell(ncell)%efield = 0.0  
            cell(ncell)%e_pot  = 0.0  
            cell(ncell)%plasma(1) = 10.0     ! Te (eV)
            cell(ncell)%plasma(2) = 10.0     ! Ti (eV)
            cell(ncell)%plasma(3) = 1.1E+12  ! ni (eV) (ne=ni assumed at present)
          ENDIF
        ENDDO

        ntri = 0       
        DEALLOCATE(be)
c        DEALLOCATE(wall_z)
c        DEALLOCATE(wall_ik)
c        DEALLOCATE(wall_ir)

      ELSE
c       -----------------------------------------------------------------
c...    Standard fluid grid:

c...    Rough crack at e-potential:
        e_pot = 0.0
        DO ir = irsep, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          e_pot(0,ir) = 0.0
          DO ik = 1, nks(ir)
            IF (ik.EQ.1) THEN
              e_pot(ik,ir) = -0.5 * kes(ik,ir) * kss(ik,ir)
            ELSE
              e_pot(ik,ir) = e_pot(ik-1,ir) -
     .                       kes(ik,ir) * (kss(ik,ir) - kss(ik-1,ir))
            ENDIF
            IF (ik.EQ.nks(ir)) e_pot(ik+1,ir) = e_pot(ik,ir) -
     .                           kes(ik,ir) * (ksmaxs(ir) - kss(ik,ir))
          ENDDO
c          DO ik = 0, nks(ir)+1
c            WRITE(0,*) 'E_POT:',ik,ir,e_pot(ik,ir)
c          ENDDO  
          DO ik = 1, nks(ir)
            frac = kss(ik,ir) / ksmaxs(ir)
            e_pot(ik,ir) = e_pot(ik,ir) - frac * e_pot(nks(ir)+1,ir)
          ENDDO  
          e_pot(nks(ir)+1,ir) = 0.0
c          DO ik = 0, nks(ir)+1
c            WRITE(0,*) 'E_POT:',ik,ir,e_pot(ik,ir)
c          ENDDO  
c          STOP 'sdfsd'
        ENDDO
        ir = irsep
        maxepot = -1.0E+20
        DO ik = 1, nks(ir)
          IF (e_pot(ik,ir).GT.maxepot) maxepot = e_pot(ik,ir)
        ENDDO
        DO ir = 2, irsep-1
          e_pot(1:nks(ir),ir) = maxepot
        ENDDO
        
c...    Load cell geometry and volume quantities for the magnetic/fluid grid:
c        fact = qtim * qtim * emi / crmi
        fact = 1.0
        ncell = 0
        DO ir = 1, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          IF (ir.LT.irsep) THEN
            ike = nks(ir) - 1
          ELSE
            ike = nks(ir)
          ENDIF
          DO ik = 1, ike
            ncell = ncell + 1
            cell(ncell)%index     = 0
            cell(ncell)%sideindex = 0
            cell(ncell)%type = 1
c...        Cell indices on the magnetic grid:
            cell(ncell)%index(1) = ik                               ! Cell index
            cell(ncell)%index(2) = ir                               ! Ring index
c...        Radial cell surfaces:
            cell(ncell)%sideindex(1,2) = 23                         ! Side index for cell surface 2
            cell(ncell)%sideindex(1,4) = 14                         ! Side index for cell surface 4
c...        Poloidal surfaces of note (targets):
            IF (ir.GE.irsep) THEN
              IF     (ik.EQ.1      ) THEN
                cell(ncell)%sideindex(2,1) = IKLO        ! Target index for cell surface 1
                cell(ncell)%sideindex(5,1) = nimindex(idds(ir,2))  ! Map to DIVIMP wall array
              ELSEIF (ik.EQ.nks(ir)) THEN
                cell(ncell)%sideindex(2,3) = IKHI  
                cell(ncell)%sideindex(5,3) = nimindex(idds(ir,1))
              ENDIF
            ENDIF
c...        Cell vertices (4 vertices assumed at present):            
            id = korpg(ik,ir)
            DO i1 = 1, nvertp(id)
              IF (ALLOCATED(d_rvertp)) THEN
                cell(ncell)%r(i1) = d_rvertp(i1,id)       ! x (m)
                cell(ncell)%z(i1) = d_zvertp(i1,id)       ! y
              ELSE
                cell(ncell)%r(i1) = DBLE(rvertp(i1,id))   ! x (m)
                cell(ncell)%z(i1) = DBLE(zvertp(i1,id))   ! y
              ENDIF
            ENDDO
c...        B-field components (approximate):         
            x(1) = DBLE(0.5 * (rvertp(1,id) + rvertp(2,id)))        ! x midpoint of the poloidal cell surface 1
            y(1) = DBLE(0.5 * (zvertp(1,id) + zvertp(2,id)))        ! y midpoint
            x(2) = DBLE(0.5 * (rvertp(3,id) + rvertp(4,id)))        ! x midpoint of the poloidal cell surface 3
            y(2) = DBLE(0.5 * (zvertp(3,id) + zvertp(4,id)))        ! y midpoint
            deltax = (x(2) - x(1))
            deltay = (y(2) - y(1))
            brat = DBLE(bratio(ik,ir))                              ! B_poloidal / B_total
            IF (DABS(deltay).LT.1.0D-10) THEN
              beta = 0.0D0
            ELSE
              beta = deltax / deltay
            ENDIF 
            Bz = DSQRT(1.0 - brat**2)                               
            By = brat * DSQRT(1.0D0/(1.0D0+beta**2))*DSIGN(1.0D0,deltay)
            Bx = beta * By
c...        CBPHI is the on-axis B-field value specified in the OSM input file:
            Bfrac = cbphi * r0 / rs(ik,ir)                         ! Rough scaling of B-field
c...        Plasma quantities:
            cell(ncell)%plasma(1) = ktebs(ik,ir)                   ! Te (eV)
            cell(ncell)%plasma(2) = ktibs(ik,ir)                   ! Ti (eV)
            cell(ncell)%plasma(3) = knbs (ik,ir)                   ! ni (eV) (ne=ni assumed at present)
            cell(ncell)%plasma(4) = SNGL(Bx) * kvhs(ik,ir)         ! vx (m-1 s-1)
            cell(ncell)%plasma(5) = SNGL(By) * kvhs(ik,ir)         ! vy
            cell(ncell)%plasma(6) = SNGL(Bz) * kvhs(ik,ir)         ! vz
c...        E&M quantities:
            cell(ncell)%bfield(1) = SNGL(Bx) * Bfrac               ! Bx (Tesla) (normalized on Eirene side)
            cell(ncell)%bfield(2) = SNGL(By) * Bfrac               ! By 
            cell(ncell)%bfield(3) = SNGL(Bz) * Bfrac               ! Bz 
            cell(ncell)%bfield(4) = brat                           ! Bratio 
            cell(ncell)%efield(1) = SNGL(Bx) * kes(ik,ir) / fact   ! Ex (not required by EIRENE)
            cell(ncell)%efield(2) = SNGL(By) * kes(ik,ir) / fact   ! Ey
            cell(ncell)%efield(3) = SNGL(Bz) * kes(ik,ir) / fact   ! Ez
        
            cell(ncell)%e_pot = e_pot(ik,ir)                       ! Electric potential (estimate)
          ENDDO
        ENDDO
c       -----------------------------------------------------------------
      ENDIF

c...  Special option to reduce the plasma density to a very low level when
c     looking at gas dynamics:
      IF (whipe.GT.0) cell(1:ncell)%plasma(3) = 1.0E+12

      IF (tetrahedrons) THEN
c   *** HACK ***
c        cell(1)%plasma(1) = 10.0                  ! Te (eV)
c        cell(1)%plasma(2) = 10.0                  ! Ti (eV)
c        cell(1)%plasma(3) = 5.0E+21               ! ni (eV) (ne=ni assumed at present)
        fact = ECH
      ELSE
        fact = ECH 
      ENDIF
      IF (eirsrcmul.NE.1.0) STOP 'STOP: IS EIRSRCMUL CODED CORRECTLY?'

c...  Specify target plasma quantities:
      it = 0
      DO region = IKLO, IKHI
        DO ir = irsep, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          IF (region.EQ.IKLO) THEN          
            ik = 1
            in = idds(ir,2)
          ELSE
            ik = nks(ir)
            in = idds(ir,1)
          ENDIF

          it = it + 1
          tardat(it,1)  = REAL(region)                                    ! Target index (IKLO=1)
          tardat(it,2)  = REAL(ik)                                        ! Cell index on fluid grid (EDGE2D: not sure how this works...)
          tardat(it,3)  = REAL(ir)                                        ! Ring index
          tardat(it,4)  = REAL(0)                                         ! Dummy, for 3D when it comes...
          tardat(it,5)  = knds(in)                                        ! ne (m-3) (required)
          tardat(it,6)  = kteds(in)                                       ! Te (eV)  (required)
          tardat(it,7)  = GetFlux(region,ir) * fact                       ! ion flux to target for species 1 (Amps) (req.)
          tardat(it,8)  = ktids(in)                                       ! Ti (eV)  (required)
          tardat(it,9)  = knds(in)                                        ! ni (m-3) (required)
          tardat(it,10) = kvds(in)                                        ! v_para (m s-1)
          tardat(it,11) = GetMach(kvds(in),kteds(in),ktids(in))           ! Mach no.
          tardat(it,12) = GetJsat(kteds(in),ktids(in),knds(in),kvds(in))  ! jsat 
          tardat(it,13) = costet(in)                                      ! *** HACK *** COSTET, for calculating tetrahedron cell particle flux...

c          WRITE(88,*) 'targ flux',ik,ir,region,GetFlux(region,ir),fact
        ENDDO
      ENDDO      
      ntardat = it

c      WRITE(0,*) 'DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: SetupEireneStrata
c
c
      SUBROUTINE SetupEireneStrata
      USE mod_eirene06
      USE mod_osm_input
      USE mod_sol28_global
      use mod_params
      use mod_slcom
      use mod_pindata
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'slcom'
c     include 'pindata'

      INTEGER i1,type,nspez,nasor,i2,insor,is,target,srcsrf
      LOGICAL assign_LO,assign_HI,assign_volrec

      TYPE(type_strata) :: tmpstrata      

      nstrata = 0
      srcsrf  = 0 

c...  Decide if default strata should be assigned:
c      assign_LO     = .TRUE.
c      assign_HI     = .TRUE.
c      assign_volrec = .TRUE.
c      DO is = 1, opt_eir%nstrata
c        IF (NINT(opt_eir%type(is))  .EQ.1.AND.
c     .           opt_eir%target(is).EQ.IKLO) assign_LO     =.FALSE.  ! -1 = LO target, -2 = high target
c        IF (NINT(opt_eir%type(is))  .EQ.1.AND.
c     .           opt_eir%target(is).EQ.IKHI) assign_HI     =.FALSE.
c        IF (     opt_eir%type(is)   .EQ.2.0) assign_volrec =.FALSE.
c      ENDDO
      assign_LO     = .FALSE.
      assign_HI     = .FALSE.
      assign_volrec = .FALSE.

c      WRITE(0,*) 'STRATA:',assign_LO,assign_HI,assign_volrec


c...  Low IK target:
      IF (assign_LO) THEN
        srcsrf  = srcsrf + 1 
        nstrata = nstrata + 1

        strata(nstrata)%range_cell(1) = -1  ! No idea what this is for, also there for ASSIGN_HI...

        strata(nstrata)%type    = 1.0
        strata(nstrata)%indsrc  = 1
        strata(nstrata)%txtsou  = '* D+ bulk ions, low index target'
c        strata(nstrata)%npts    = 100
        strata(nstrata)%npts    = -90000
c
c       jdemod - ninitl appears to be the random number seed for eirene
c              - presumably -1 means to generate a new seed?
c              - use piniseed if > 0
c
        if (piniseed.le.0) then 
           strata(nstrata)%ninitl  = -1
        else
           strata(nstrata)%ninitl  = piniseed
        endif
c
        strata(nstrata)%nemods  =  3
        strata(nstrata)%flux    = 1.0
        strata(nstrata)%species_tag = 'FFFT'
        strata(nstrata)%nspez   =  1
        strata(nstrata)%distrib = 'FFTFF'
        strata(nstrata)%inum    =  1
        strata(nstrata)%indim   =  4
        strata(nstrata)%insor   = -srcsrf  ! -2
        strata(nstrata)%sorwgt  = 1.0
        strata(nstrata)%sorlim  = 124.0
        strata(nstrata)%sorind  = 1.0
        strata(nstrata)%nrsor   = -1
        strata(nstrata)%nasor   =  0
        strata(nstrata)%sorad   = 0.0
        strata(nstrata)%sorene  = 3.0
        strata(nstrata)%soreni  = 0.5
        strata(nstrata)%sorcos  = 1.0
        strata(nstrata)%sormax  = 90.0
      ENDIF

c...  High IK target:
      IF (assign_HI) THEN
        srcsrf  = srcsrf + 1 
        nstrata = nstrata + 1

        strata(nstrata)%range_cell(1) = -2  ! Again, can't recall what I was starting to setup here, 
c        strata(nstrata)%range_cell(1) = 1  ! but I think it was something for Lorenzo and the ULS...
c        strata(nstrata)%range_cell(2) = 1
c        strata(nstrata)%range_tube(1) = 1
c        strata(nstrata)%range_tube(2) = 1E+6

        strata(nstrata)%type    = 1.0
        strata(nstrata)%indsrc  = 1
        strata(nstrata)%txtsou  = '* D+ bulk ions, high index target'
c        strata(nstrata)%npts    = 100
        strata(nstrata)%npts    = -90000
c
c       jdemod - ninitl appears to be the random number seed for eirene
c              - presumably -1 means to generate a new seed?
c              - use piniseed if > 0
c
        if (piniseed.le.0) then 
           strata(nstrata)%ninitl  = -1
        else
           strata(nstrata)%ninitl  = piniseed
        endif
c
        strata(nstrata)%nemods  =  3
        strata(nstrata)%flux    = 1.0
        strata(nstrata)%species_tag = 'FFFT'
        strata(nstrata)%nspez   =  1
        strata(nstrata)%distrib = 'FFTFF'
        strata(nstrata)%inum    =  1
        strata(nstrata)%indim   =  4
        strata(nstrata)%insor   = -srcsrf  ! -3
        strata(nstrata)%sorwgt  = 1.0
        strata(nstrata)%sorlim  = 124.0
        strata(nstrata)%sorind  = 2.0
        strata(nstrata)%nrsor   = -1
        strata(nstrata)%nasor   =  0
        strata(nstrata)%sorad   = 0.0
        strata(nstrata)%sorene  = 3.0
        strata(nstrata)%soreni  = 0.5
        strata(nstrata)%sorcos  = 1.0
        strata(nstrata)%sormax  = 90.0
      ENDIF

c...  Volume recombination:
      IF (assign_volrec) THEN
        nstrata = nstrata + 1
        strata(nstrata)%type    = 2.0
        strata(nstrata)%indsrc  = 1
        strata(nstrata)%txtsou  = '* Volume recombination'
c        strata(nstrata)%npts    = 100
        strata(nstrata)%npts    = -90000
c
c       jdemod - ninitl appears to be the random number seed for eirene
c              - presumably -1 means to generate a new seed?
c              - use piniseed if > 0
c
        if (piniseed.le.0) then 
           strata(nstrata)%ninitl  = -1
        else
           strata(nstrata)%ninitl  = piniseed
        endif
c
        strata(nstrata)%nemods  =  3
        strata(nstrata)%flux    = 1.0
        strata(nstrata)%species_tag = 'FFFT'
        strata(nstrata)%nspez   =  1
        strata(nstrata)%distrib = 'FFFTF'
        strata(nstrata)%inum    =  1
        strata(nstrata)%indim   =  4
        strata(nstrata)%insor   =  5
        strata(nstrata)%nasor   =  0
        strata(nstrata)%sorwgt  = 1.0
        strata(nstrata)%sorlim  = 124.0
        strata(nstrata)%sorind  = 0.0
        strata(nstrata)%nrsor   = -1
        strata(nstrata)%sorad   = 0.0
        strata(nstrata)%sorene  = 3.0
        strata(nstrata)%soreni  = 0.0
        strata(nstrata)%sorcos  = 1.0
        strata(nstrata)%sormax  = 90.0
      ENDIF

c...  User specified neutral injection/puffing:
      DO is = 1, opt_eir%nstrata
        WRITE(PINOUT,*) 'STRATA:TYPE=',NINT(opt_eir%type(is))

        SELECTCASE (NINT(opt_eir%type(is)))
          CASE (999)  ! For compatibility with old strata definition format, see below
          CASE (1)
            srcsrf = srcsrf + 1 
            target = opt_eir%target(is)
            nstrata = nstrata + 1
            strata(nstrata)%range_cell(1) = -1  ! No idea what this is for...
            strata(nstrata)%type    = 1.0
            strata(nstrata)%indsrc  = 1
            strata(nstrata)%txtsou  = '* '//opt_eir%txtsou(is)
            strata(nstrata)%npts    = opt_eir%npts(is)
c
c           jdemod - ninitl appears to be the random number seed for eirene
c              - presumably -1 means to generate a new seed?
c              - use piniseed if > 0
c
            if (piniseed.le.0) then 
               strata(nstrata)%ninitl  = -1
            else
               strata(nstrata)%ninitl  = piniseed
            endif
c
            strata(nstrata)%nemods  =  3
            strata(nstrata)%flux    = opt_eir%flux(is)
            strata(nstrata)%species_tag = 'FFFT'
            strata(nstrata)%nspez   =  1
            strata(nstrata)%distrib = 'FFTFF'
            strata(nstrata)%inum    =  1
            strata(nstrata)%indim   =  4
            strata(nstrata)%insor = -srcsrf
c            IF (target.EQ.IKLO) strata(nstrata)%insor = -2
c            IF (target.EQ.IKHI) strata(nstrata)%insor = -3
            strata(nstrata)%sorwgt  = 1.0
            strata(nstrata)%sorlim  = 124.0
            strata(nstrata)%sorind  = REAL(srcsrf)
c            IF (target.EQ.IKLO) strata(nstrata)%sorind = 1.0
c            IF (target.EQ.IKHI) strata(nstrata)%sorind = 2.0
            strata(nstrata)%nrsor   = -1
            strata(nstrata)%nasor   =  0
            strata(nstrata)%sorad   = 0.0
            strata(nstrata)%sorene  = 3.0
            strata(nstrata)%soreni  = 0.5
            strata(nstrata)%sorcos  = 1.0
            strata(nstrata)%sormax  = 90.0
          CASE (2)
            nstrata = nstrata + 1
            strata(nstrata)%type    = 2.0
            strata(nstrata)%indsrc  = 1
            strata(nstrata)%txtsou  = '* '//opt_eir%txtsou(is)
            strata(nstrata)%npts    = opt_eir%npts(is)
c
c           jdemod - ninitl appears to be the random number seed for eirene
c              - presumably -1 means to generate a new seed?
c              - use piniseed if > 0
c
            if (piniseed.le.0) then 
               strata(nstrata)%ninitl  = -1
            else
               strata(nstrata)%ninitl  = piniseed
            endif
c
            strata(nstrata)%nemods  =  3
            strata(nstrata)%flux    = opt_eir%flux(is)
            strata(nstrata)%species_tag = 'FFFT'
            strata(nstrata)%nspez   =  1
            strata(nstrata)%distrib = 'FFFTF'
            strata(nstrata)%inum    =  1
            strata(nstrata)%indim   =  4
            strata(nstrata)%insor   =  5
            strata(nstrata)%nasor   =  0
            strata(nstrata)%sorwgt  = 1.0
            strata(nstrata)%sorlim  = 124.0
            strata(nstrata)%sorind  = 0.0
            strata(nstrata)%nrsor   = -1
            strata(nstrata)%sorad   = 0.0
            strata(nstrata)%sorene  = 3.0
            strata(nstrata)%soreni  = 0.0
            strata(nstrata)%sorcos  = 1.0
            strata(nstrata)%sormax  = 90.0
          CASE (3)

            write(0,*) 'puff index weirdness',is,nstrata+1, 
     .                 '  '//TRIM(opt_eir%txtsou(is))

            IF (opt_eir%type(is).EQ.3.1) THEN  ! *** (small) HACK ***
              nspez = 2
              insor = 2 ! 1
              nasor = 1
              beam  = 1 ! Turn on beams in EIRENE input file...
            ELSE
              nspez = 1
              insor = 1
              nasor = 0
            ENDIF
            nstrata = nstrata + 1
            strata(nstrata)%type    = opt_eir%type(is)
            strata(nstrata)%indsrc  = -1
            strata(nstrata)%txtsou  = '* point injection, '//
     .                                opt_eir%txtsou(is)
c
c           jdemod - ninitl appears to be the random number seed for eirene
c              - presumably -1 means to generate a new seed?
c              - use piniseed if > 0
c
            if (piniseed.le.0) then 
               strata(nstrata)%ninitl  = -1
            else
               strata(nstrata)%ninitl  = piniseed
            endif
c
            strata(nstrata)%nemods  =  1
            strata(nstrata)%npts    = opt_eir%npts(is)
            strata(nstrata)%flux    = opt_eir%flux(is) *
     .                                opt_eir%flux_fraction(is)
            strata(nstrata)%species_tag = 'FFFF'
            i2 = opt_eir%species(is)
            strata(nstrata)%species_tag(i2:i2) = 'T'
            strata(nstrata)%nspez   = nspez
            strata(nstrata)%distrib = 'TFFFF'
            strata(nstrata)%inum    = 0
            strata(nstrata)%indim   = 0
            strata(nstrata)%insor   = insor
            strata(nstrata)%sorwgt  = 1.0
            strata(nstrata)%sorlim  = 220.0
            strata(nstrata)%sorind  = 0.0
            strata(nstrata)%nrsor   = 0
            strata(nstrata)%nasor   = nasor
            strata(nstrata)%sorad(1:6) =opt_eir%sorad(1:6,is) *  ! Convert to cm
     .                                  100.0  
c            WRITE(0,*) strata(nstrata)%sorad(1:3)
c            STOP 'sdgsdgd'
            strata(nstrata)%sorene  = opt_eir%sorene(is)
            strata(nstrata)%soreni  = opt_eir%sorene(is)  ! bug - seems the puff energy is set from SORENI, not SORENE, so all puff had E=0 before this? -SL, 19/01/12
            strata(nstrata)%sorcos  = opt_eir%sorcos(is)
            strata(nstrata)%sormax  = opt_eir%sormax(is)
          CASE (4)  ! Puff wall surface (not working...)
            nspez = 1
            insor = 1
            nasor = 0

            nstrata = nstrata + 1
            strata(nstrata)%type    = opt_eir%type(is)
            strata(nstrata)%indsrc  = -1
            strata(nstrata)%txtsou  = '* surface injection, '//
     .                                opt_eir%txtsou(is)
c
c           jdemod - ninitl appears to be the random number seed for eirene
c              - presumably -1 means to generate a new seed?
c              - use piniseed if > 0
c
            if (piniseed.le.0) then 
               strata(nstrata)%ninitl  = -1
            else
               strata(nstrata)%ninitl  = piniseed
            endif
c
            strata(nstrata)%nemods  =  1
            strata(nstrata)%npts    = opt_eir%npts(is)
            strata(nstrata)%flux    = opt_eir%flux(is) *
     .                                opt_eir%flux_fraction(is)
            strata(nstrata)%species_tag = 'FFFF'
            i2 = opt_eir%species(is)
            strata(nstrata)%species_tag(i2:i2) = 'T'
            strata(nstrata)%nspez   = nspez
            strata(nstrata)%distrib = 'FFTFF'
            strata(nstrata)%inum    = 0
            strata(nstrata)%indim   = 1 ! 0
            strata(nstrata)%insor   = 2 ! insor  - surface in 3a list
            strata(nstrata)%sorwgt  = 1.0
            strata(nstrata)%sorlim  = 202.0  ! 220.0
            strata(nstrata)%sorind  = 0.0
            strata(nstrata)%nrsor   = 0
            strata(nstrata)%nasor   = nasor
            strata(nstrata)%sorad(1:6) =opt_eir%sorad(1:6,is) *  ! Convert to cm
     .                                  100.0  
            strata(nstrata)%sorene  = opt_eir%sorene(is)
            strata(nstrata)%soreni  = opt_eir%sorene(is) 

            write(0,*) 'puff eni',opt_eir%sorene(is)
            strata(nstrata)%sorcos  = opt_eir%sorcos(is)
            strata(nstrata)%sormax  = opt_eir%sormax(is)
c            WRITE(0,*) 'WALL SURFACE NEUTRAL INJECTION NOT WORKING'
c            STOP 
          CASE DEFAULT
            CALL ER('SetupEireneStrata','Unrecognized strata type',*99)
        ENDSELECT

      ENDDO

c...  OLD SPECIFICATION: User specified neutral injection/puffing:
      DO i1 = 1, eirnpuff
        STOP 'OLD PUFF SPECIFICATION OBSOLETE'
        SELECTCASE (NINT(eirpuff(1,i1)))
          CASE (999)
          CASE (1)
            IF (eirpuff(1,i1).EQ.1.1) THEN  !***(small) hack***
              nspez = 2
              insor = 2 ! 1
              nasor = 1
              beam  = 1 ! Turn on beams in input file...
            ELSE
              nspez = 1
              insor = 1
              nasor = 0
            ENDIF
            eirpuff(1,i1) = 999.0
            nstrata = nstrata + 1
            strata(nstrata)%indsrc  = -1
            strata(nstrata)%type    = 999.0
            strata(nstrata)%txtsou  = '* point injection, '//
     .                                eircpuff(i1)
            strata(nstrata)%npts    = NINT(eirpuff(2,i1))
c
c           jdemod - ninitl appears to be the random number seed for eirene
c              - presumably -1 means to generate a new seed?
c              - use piniseed if > 0
c
            if (piniseed.le.0) then 
               strata(nstrata)%ninitl  = -1
            else
               strata(nstrata)%ninitl  = piniseed
            endif
c
            strata(nstrata)%nemods  =  1
            strata(nstrata)%flux    = eirpuff(3,i1) * eirpuff(4,i1)
            strata(nstrata)%species_tag = 'FFFF'
            i2 = NINT(eirpuff(5,i1))
            strata(nstrata)%species_tag(i2:i2) = 'T'
            strata(nstrata)%nspez   =  nspez
            strata(nstrata)%distrib = 'TFFFF'
            strata(nstrata)%inum    =  0
            strata(nstrata)%indim   =  0
            strata(nstrata)%insor   =  insor
            strata(nstrata)%sorwgt  = 1.0
            strata(nstrata)%sorlim  = 220.0
            strata(nstrata)%sorind  = 0.0
            strata(nstrata)%nrsor   =  0
            strata(nstrata)%nasor   =  nasor
            strata(nstrata)%sorad(1:6) = eirpuff(9:14,i1)
            strata(nstrata)%sorad(1:3) = strata(nstrata)%sorad(1:3) *  ! Convert to cm
     .                                   100.0  
            strata(nstrata)%sorene  = eirpuff(6,i1)
            strata(nstrata)%soreni  = 0.0
            strata(nstrata)%sorcos  = eirpuff(7,i1)
            strata(nstrata)%sormax  = eirpuff(8,i1)
          CASE (2)  ! Puff from non-standard default surface (not working...)
            nstrata = nstrata + 1
            strata(nstrata)%indsrc  = -1
            strata(nstrata)%type    = -999.0
            strata(nstrata)%txtsou  = '* non-default surface puff, '//
     .                                eircpuff(i1)
            strata(nstrata)%npts    = NINT(eirpuff(2,i1))
c
c           jdemod - ninitl appears to be the random number seed for eirene
c              - presumably -1 means to generate a new seed?
c              - use piniseed if > 0
c
            if (piniseed.le.0) then 
               strata(nstrata)%ninitl  = -1
            else
               strata(nstrata)%ninitl  = piniseed
            endif
c
            strata(nstrata)%nemods  =  1
            strata(nstrata)%flux    = eirpuff(3,i1) * eirpuff(4,i1)
            strata(nstrata)%species_tag = 'FFFF'
            i2 = NINT(eirpuff(5,i1))
            strata(nstrata)%species_tag(i2:i2) = 'T'
            strata(nstrata)%nspez   = 1
            strata(nstrata)%distrib = 'FFTFF'
            strata(nstrata)%inum    = 0
            strata(nstrata)%indim   = NINT(eirpuff(9 ,i1))
            strata(nstrata)%insor   = NINT(eirpuff(10,i1))
            strata(nstrata)%sorwgt  = 1.0
            strata(nstrata)%sorlim  = eirpuff(11,i1)
            strata(nstrata)%sorind  = 0.0
            strata(nstrata)%nrsor   =  0
            strata(nstrata)%nasor   = NINT(eirpuff(12,i1)) 
c            strata(nstrata)%sorad =  0.0
            strata(nstrata)%sorad(1) =  1.94E+01
            strata(nstrata)%sorad(2) =  1.94E+01 
            strata(nstrata)%sorad(3) = -100.0
            strata(nstrata)%sorad(4) =  100.0
            strata(nstrata)%sorad(5) =  0.0
            strata(nstrata)%sorad(6) =  0.0
            strata(nstrata)%sorene  = eirpuff(6,i1)
            strata(nstrata)%soreni  = 0.0
            strata(nstrata)%sorcos  = eirpuff(7,i1)
            strata(nstrata)%sormax  = eirpuff(8,i1)
          CASE DEFAULT
        ENDSELECT
      ENDDO


c...  TEMP: make sure regular strata come first until surface to strata mapping is done properly...
c      DO is = 1, nstrata
c        IF (strata(is)%type.GE.3.0) THEN
c          tmpstrata = strata(is)
c          DO i1 = is, nstrata-1
c            strata(i1) = strata(i1+1)
c          ENDDO
c          strata(nstrata) = tmpstrata  
c        ENDIF
c      ENDDO

      WRITE(PINOUT,*) 'NSTRATA:',nstrata,strata(1:nstrata)%type

c      WRITE(0,*) 'STRATA:',nstrata,strata(1:nstrata)%type
c      WRITE(0,*) 'STRATA:',strata(1:nstrata)%type
c      WRITE(0,*) 'STRATA:',alloc

c      STOP 'sdgsgsdsdgsgd'

c...  Debugging:
c      IF (time_dependent) THEN
c        strata(1:3)%npts = 10
c      ENDIF
c      strata(1:3)%npts = 100000
c      strata(1:3)%npts = 2
c      strata(1:10)%ninitl = 11111 ! 22222 ! 99887


c...  Count the number of gas strata
      gas_n = 0
      DO is = 1, nstrata
        IF (strata(is)%type.GE.3.0) gas_n = gas_n + 1
        IF (strata(is)%type.NE.3.0.AND.gas_n.GT.0) THEN
          CALL ER('SetupEireneStrata','Very sorry, but all non-gas '//
     .            'puff strata need to appear before the gas puff '//
     .            'strata in the list; see NEUTRAL SOURCES in the '//
     .            'input file',*99)
        ENDIF
      ENDDO
      WRITE(0,*) 'GAS STRATA COUNTED=',gas_n


      RETURN
 99   WRITE(0,*) '  INDEX =',is
      WRITE(0,*) '  TYPE  =',strata(is)%type
      STOP
      END
c
c ======================================================================
c ======================================================================
c
      SUBROUTINE ReadEireneResults_06(iitersol)
      USE mod_eirene06
      USE mod_sol28_global
      IMPLICIT none

      INTEGER iitersol

      CALL ReadParticleTracks_04

      CALL LoadEireneData_06(iitersol,opt_eir%ilspt)
     
      call targflux
c
c     Print out some of the Eirene data in DIVIMP format
c
      call eireprn

      CALL WriteEireneDatFileInfo

c...  Move this back to near where the arrays are allocated, so that the memory is
c     not occupied when EIRENE is actually running, although will need to get
c     the cell index mapping back from EIRENE, which isn't happening at present...
      CALL DEALLOC_VERTEX  
      CALL DEALLOC_SURFACE   
      CALL DEALLOC_TRIANGLE

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
c
      SUBROUTINE LoadEireneData_06(iitersol,ilspt)
      USE mod_interface
      USE mod_geometry
      USE mod_eirene06_parameters
      USE mod_eirene06
      USE mod_eirene06_locals
      USE mod_eirene_history
      USE mod_divimp
      USE mod_options
      use mod_params
      use mod_comtor
      use mod_cgeom
      use mod_pindata
      use mod_slcom
      IMPLICIT none
 
      INTEGER, INTENT(IN) :: iitersol,ilspt

c     INCLUDE 'params'
c     INCLUDE 'comtor'
c     INCLUDE 'cgeom'
c     INCLUDE 'pindata'
c     INCLUDE 'slcom'

      INTEGER fp,ntally,ndata,icount,index(30),ik,ir,i1,i2,iside,isrf,
     .        iblk,iatm,imol,iion,ipho,ilin,isur,cvesm(MAXSEG),tube,ike,
     .        i3,i4,idum,i,in,iw,j,ivtx,iliin,code
      LOGICAL goodeof,debug,wall_ignored,warning_message,binary
      REAL    rdum(30),frac,norm,len,cir,area,volume,fact,
     .        sumion,amps,pflux,net1,net2,v1(3),v2(3)
      CHARACTER buffer*256,species*32,fname*1024,tag*64,snip*6

      INTEGER iobj,igrp
      INTEGER*2, ALLOCATABLE :: fluid_ik(:),fluid_ir(:),wall_in(:,:)
      INTEGER  , ALLOCATABLE :: vmap(:,:)
      REAL     , ALLOCATABLE :: tdata(:,:,:),tflux(:,:),eirdat(:,:,:),
     .                          sflux(:,:)

      write(0,*) 'PIN:'//TRIM(actpin)//':'

      IF (TRIM(actpin).EQ.' reire06') THEN
        WRITE(0,*) 'BINARY EIRENE DATA READ ACTIVATED'
        binary = .TRUE.
      ELSE
        binary = .FALSE.
      ENDIF

      wall_ignored = .FALSE.
      debug        = .FALSE.
      icount = 0
      iside  = 0

      ALLOCATE(tdata(MAXNKS,MAXNRS,5 ))
      ALLOCATE(tflux(MAXSEG       ,14))
      tdata = 0.0
      tflux = 0.0

      IF (ilspt.GT.0) THEN
        ALLOCATE(eirdat(MAXNKS,MAXNRS,2))
        eirdat = 0.0
      ENDIF



      net1 = 0.0
      net2 = 0.0




      IF (tetrahedrons) THEN
        ALLOCATE(fluid_ik(nobj))
        ALLOCATE(fluid_ir(nobj))
        ALLOCATE(wall_in(4,nobj))
        ALLOCATE(sflux(nsrf,4))
        fluid_ik = 0
        fluid_ir = 0
        sflux = 0.0
        DO iobj = 1, nobj
          igrp = obj(iobj)%group
          IF (grp(igrp)%origin.EQ.GRP_MAGNETIC_GRID) THEN
            fluid_ik(iobj) = obj(iobj)%index(IND_IK)
            fluid_ir(iobj) = obj(iobj)%index(IND_IR)
          ENDIF
          DO iside = 1, 4
            isrf = ABS(obj(iobj)%iside(iside))
            wall_in(iside,iobj) = srf(isrf)%index(IND_WALL_STD)
          ENDDO
        ENDDO
      ELSE
        ALLOCATE(fluid_ik(  ntri))
        ALLOCATE(fluid_ir(  ntri))
        ALLOCATE(wall_in (3,ntri))
        fluid_ik = 0
        fluid_ir = 0
        wall_in  = 0
        DO iobj = 1, ntri
          IF     (tri(iobj)%type.EQ.MAGNETIC_GRID) THEN  
            fluid_ik(iobj) = tri(iobj)%index(1)                   ! Should pull these from .transfer
            fluid_ir(iobj) = tri(iobj)%index(2)   
            ik = fluid_ik(iobj)
            ir = fluid_ir(iobj)
c            IF (ir.LT.irsep) CYCLE
c            IF (ik.EQ.1      ) wall_in(1:3,iobj) = nimindex(idds(ir,2))  ! A problem in corners..?
c            IF (ik.EQ.nks(ir)) wall_in(1:3,iobj) = nimindex(idds(ir,1))
          ELSEIF (tri(iobj)%type.EQ.VACUUM_GRID) THEN
c            wall_in(1:3,iobj) = tri(iobj)%sideindex(3,1:3)  ! =surface(i2)%index(1) 
c            IF (iobj.EQ.15033) THEN
c              WRITE(0,*) 'wall_in',wall_in(:,iobj)
c              WRITE(0,*) 'wall_in',tri(iobj)%sideindex(3,:)
c            ENDIF
          ENDIF
          wall_in(1:3,iobj) = tri(iobj)%sideindex(5,1:3)
        ENDDO
      ENDIF

      goodeof = .FALSE.

c     - Impurity momentum source [(kg x m)/(cell x sec^2)]
c     - H neutral density        [particles/cell]             PINATOM
c     - Impurity neutral density [particles/cell]             PINZ0
c     - H2 density               [molecules/cell]             PINMOL
c...Changed on Dec 10, 1999 - SL
c     - H2+ density              [molecules/cell]             PINMOI
c OLD     - H2+ density              [molecules/cell]            +PINMOL
c     - Avg. neutr. H energy     [eV]                         PINENA
c     - Avg. neutr. imp. energy  [ev]                         PINENZ
c     - Avg. H2 energy           [eV]                         PINENM
c     - Avg. H2+ energy          [eV]  !!! (not yet added to DIVIMP) !!!
c     - Halpha emissivity        [photons/(cell x sec)]       PINALPHA
c       (5 blocks for different generation mechanisms)
c     - H rec. source            [particles/(cell x sec)]     PINREC
c     - H+ (+imp.) energy source [Watt/cell]                  PINQI
c     - Electron energy source   [Watt/cell]                  PINQE

      pinalpha = 0.0
      pinline  = 0.0
      pinatom  = 0.0
      pinena   = 0.0
      pinmol   = 0.0
      pinenm   = 0.0
      pinenz   = 0.0

      hescpd = 0.0
      hescal = 0.0

      IF (rel_opt.EQ.1.OR.rel_opt.EQ.3) THEN
        frac = rel_frac
      ELSE
        frac = 1.0
      ENDIF
      IF (sldebug.NE.0)
     .  WRITE(0     ,*) 'RELAXATION FRACTION FOR EIRENE06:',frac
      WRITE(PINOUT,*) 'RELAXATION FRACTION FOR EIRENE06:',frac

c...  Open the EIRENE data transfer file:

      fp = 99
      IF (citersol.GT.0.AND.
     .    (tetrahedrons.OR.time_dependent)) THEN
        WRITE(fname,'(A,I3.3,A)') 'eirene.',iitersol,'.transfer'
      ELSE
        fname='eirene.transfer'
      ENDIF
      IF (binary) THEN
        OPEN(UNIT=fp,FILE=TRIM(fname),FORM='UNFORMATTED',
     .       STATUS='OLD',ERR=98)
      ELSE
        OPEN(UNIT=fp,FILE=TRIM(fname),ACCESS='SEQUENTIAL',
     .       STATUS='OLD',ERR=98)
      ENDIF

c...  Read through the file:
      iblk = 0
      iatm = 0
      imol = 0
      iion = 0
      ipho = 0
      ilin = 0
      cvesm = 0
      DO WHILE (.TRUE.)
        IF (binary) THEN
          READ(fp,END=10) code
c          IF (debug) WRITE(0,*) '===CODE:',code
        ELSE
          READ(fp,'(A256)',END=10) buffer
        ENDIF
c       ----------------------------------------------------------------
        IF     ((binary.AND.code.EQ.-9999911).OR.(.NOT.binary.AND.
     .           buffer(1:22).EQ.'* BULK PARTICLES - VOL')) THEN
          IF (debug) WRITE(0,*) '===BULK PARTICLES: VOLUME TALLIES==='
          iblk = iblk + 1
          IF (binary) THEN
            READ(fp,  ERR=97) ntally
            READ(fp,  ERR=97) ndata                        
            READ(fp,  ERR=97) (index(i1),i1=1,ntally)          
          ELSE
            READ(fp,*,ERR=97) ntally
            READ(fp,*,ERR=97) ndata                         ! Check...
            READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          ENDIF
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum,binary)
            IF (fluid_ik(icount).NE.0) THEN
              ik = fluid_ik(icount)                          ! Should pull these from .transfer
              ir = fluid_ir(icount)
              IF     (iblk.EQ.1) THEN                           ! Data for D+ only:
                tdata(ik,ir,1) = tdata(ik,ir,1) + rdum(5)   ! PINION
                tdata(ik,ir,2) = tdata(ik,ir,2) + rdum(6)   ! PINREC (relax?)
                tdata(ik,ir,3) = tdata(ik,ir,3) + rdum(7)   ! PINMP
                tdata(ik,ir,4) = tdata(ik,ir,4) + rdum(8)   ! PINQi
                tdata(ik,ir,5) = tdata(ik,ir,5) + rdum(9)   ! PINQe
              ELSEIF (iblk.EQ.2.AND.ilspt.GT.0) THEN                     
                eirdat(ik,ir,1) = eirdat(ik,ir,1) + rdum(5)  ! Impurity atom (probably) ionisation
              ENDIF
            ENDIF
          ENDDO
          IF (debug) WRITE(0,*) '===DONE==='
        ELSEIF ((binary.AND.code.EQ.-9999912).OR.(.NOT.binary.AND.
     .          buffer(1:22).EQ.'* BULK PARTICLES - SUR')) THEN
          IF (debug) WRITE(0,*) '===BULK PARTICLES: SURFACE===',iblk
          IF (binary) THEN
            READ(fp  ,ERR=97) ntally
            READ(fp  ,ERR=97) ndata                         
            READ(fp  ,ERR=97) (index(i1),i1=1,ntally)          
          ELSE
            READ(fp,*,ERR=97) ntally
            READ(fp,*,ERR=97) ndata                         
            READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          ENDIF
          icount = 1
          DO WHILE (icount.LE.ndata)
            icount = icount + 1
            CALL NextLine(fp,ntally,iobj,rdum,binary)
            iside = NINT(rdum(1))
            in = wall_in(iside,iobj)
            IF (in.EQ.0) THEN
c             Lots of data comes through that's not associated with 
c             the standard DIVIMP neutral wall, ignore for now...              
c              IF (tri(iobj)%sideindex(4,iside).NE.0) wall_ignored=.TRUE.
              wall_ignored=.TRUE.
              CYCLE
            ENDIF
            IF     (iblk.EQ.1) THEN                ! Only for D, presumably the 1st atom species, need check...
              tflux(in,3) = tflux(in,3) + rdum(2)  ! Bulk particle flux (s-1)
              tflux(in,4) = tflux(in,4) + rdum(3)  ! Bulk energy flux   (eV s-1)
              cvesm(in) = 1
            ELSEIF (iblk.EQ.2.AND.ilspt.GT.0) THEN
c...          Sputtering turned on in EIRENE, ignore the data:
            ELSEIF ((iblk.EQ.2.OR.iblk.EQ.3).AND.bgk.GT.0) THEN
c...          Ignore D and D2 BGK data:
            ELSE
              CALL ER('LoadEireneData_06','IBLK out of bounds, '//
     .                'unexpected this is...',*99)
            ENDIF
c            WRITE(eirfp,*) 'STORING DATA',in,iobj,iside
          ENDDO
          IF (debug) WRITE(0,*) '===DONE==='
c       ----------------------------------------------------------------
        ELSEIF ((binary.AND.code.EQ.-9999921).OR.(.NOT.binary.AND.
     .          buffer(1:18).EQ.'* TEST ATOMS - VOL')) THEN
          iatm = iatm + 1
          IF (debug) WRITE(0,*) '===TEST ATOMS: VOLUME TALLIES===',iatm
          IF (binary) THEN
            READ(fp  ,ERR=97) ntally
            READ(fp  ,ERR=97) ndata                         
            READ(fp  ,ERR=97) (index(i1),i1=1,ntally)          
          ELSE
            READ(fp,*,ERR=97) ntally
            READ(fp,*,ERR=97) ndata                         
            READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          ENDIF
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum,binary)
            IF (fluid_ik(icount).NE.0) THEN
              ik = fluid_ik(icount)                          ! Should pull these from .transfer
              ir = fluid_ir(icount)
              IF     (iatm.EQ.1) THEN                        ! Only for D, presumably the 1st atom species, need check...
                pinatom(ik,ir) = pinatom(ik,ir) + rdum(1) 
                pinena (ik,ir) = pinena (ik,ir) + rdum(6) 
              ELSEIF (iatm.EQ.2.AND.ilspt.GT.0) THEN                     
                eirdat(ik,ir,2) = eirdat(ik,ir,2) + rdum(1)  ! Impurity atom density
                pinenz(ik,ir  ) = pinenz(ik,ir  ) + rdum(6)  ! Impurity atom energy
              ENDIF
            ENDIF
          ENDDO
          IF (debug) WRITE(0,*) '===DONE==='
c       ----------------------------------------------------------------
        ELSEIF ((binary.AND.code.EQ.-9999922).OR.(.NOT.binary.AND.
     .          buffer(1:18).EQ.'* TEST ATOMS - SUR')) THEN
          IF (debug) WRITE(0,*) '===TEST ATOMS: SURFACE FLUXES===',iatm
          IF (binary) THEN
            READ(fp  ,ERR=97) ntally
            READ(fp  ,ERR=97) ndata                         
            READ(fp  ,ERR=97) (index(i1),i1=1,ntally)          
          ELSE
            READ(fp,*,ERR=97) ntally
            READ(fp,*,ERR=97) ndata                         
            READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          ENDIF
          icount = 1
          DO WHILE (icount.LE.ndata)
            icount = icount + 1
            CALL NextLine(fp,ntally,iobj,rdum,binary)
            iside = NINT(rdum(1))
            in = wall_in(iside,iobj)
            IF (in.EQ.0) THEN
c             Lots of data comes through that's not associated with 
c             the standard DIVIMP neutral wall, ignore for now...              
c              IF (tri(iobj)%sideindex(4,iside).NE.0) wall_ignored=.TRUE.
              wall_ignored=.TRUE.
              CYCLE
            ENDIF

            IF (tetrahedrons) THEN 
              isrf = obj(iobj)%iside(iside)
              fact = SIGN(1.0,REAL(isrf))
            ELSE
              fact = 1.0
            ENDIF

            IF     (iatm.EQ.1) THEN                         ! Only for D, presumably the 1st atom species, need check...
c              WRITE(0,*) 'iobj,iside',iobj,iside
              tflux(in,1 ) = tflux(in,1 ) + fact * rdum(2)         ! Incident atom particle flux (s-1)
              tflux(in,2 ) = tflux(in,2 ) + fact * rdum(3)         ! Incident atom energy flux   (eV s-1)
              tflux(in,14) = tflux(in,14) + fact * SUM(rdum(4:7))  ! Emitted  atom particle flux (s-1)
              cvesm(in) = 1
              IF (tetrahedrons) THEN
                isrf = ABS(isrf)
                sflux(isrf,1) = 1.0
                sflux(isrf,2) = sflux(isrf,2) +
     .                          fact * (rdum(2) - SUM(rdum(4:7)))
                IF (in.EQ.13) THEN
                  net1 = net1 +  fact * (rdum(2) - SUM(rdum(4:7)))

c                  write(0,'(A,1P,2E12.4,0P,I6,1P,2E12.4,0P)') 
c     .              'net 13 atm:',net1,fact,isrf,rdum(2),SUM(rdum(4:7))
                ENDIF
              ENDIF
            ELSEIF (iatm.EQ.2.AND.ilspt.GT.0) THEN
c...          Sputtering turned on in EIRENE, assume (for now) this data 
c             is for the impurity species (loose...):
              tflux(in,5 ) = tflux(in,5 ) + rdum(4) - rdum(2)  ! Emitted atom particle flux from incident D atoms   (   s-1)
              tflux(in,6 ) = tflux(in,6 ) + rdum(8)            ! Emitted atom energy   flux          "              (eV s-1)
              tflux(in,9 ) = tflux(in,9 ) + rdum(7)            ! Emitted atom particle flux from indicent bulk ions (   s-1)
              tflux(in,10) = tflux(in,10) + rdum(9)            ! Emitted atom energy   flux           "             (eV s-1)
              tflux(in,11) = tflux(in,11) + rdum(2)            ! Incident (impurity) atom particle flux (   s-1)
              tflux(in,12) = tflux(in,12) + rdum(3)            ! Incident (impurity) atom energy   flux (eV s-1)
              IF (tetrahedrons) THEN
                isrf = ABS(isrf)
                sflux(isrf,1) = 1.0
                sflux(isrf,4) = sflux(isrf,4) +
     .                          fact * (rdum(4) - rdum(2) + rdum(7))
              ENDIF
            ELSE
              CALL ER('LoadEireneData_06','IATM out of bounds, '//
     .                'unexpected this is...',*99)
            ENDIF
c            WRITE(eirfp,*) 'STORING DATA',in,iobj,iside
          ENDDO
          IF (debug) WRITE(0,*) '===DONE==='
c       ----------------------------------------------------------------
        ELSEIF ((binary.AND.code.EQ.-9999931).OR.(.NOT.binary.AND.
     .          buffer(1:22).EQ.'* TEST MOLECULES - VOL')) THEN
          IF (debug) WRITE(0,*) '===TEST MOLECULES: VOLUME TALLIES==='
          imol = imol + 1
          IF (binary) THEN
            READ(fp  ,ERR=97) ntally
            READ(fp  ,ERR=97) ndata                         
            READ(fp  ,ERR=97) (index(i1),i1=1,ntally)          
          ELSE
            READ(fp,*,ERR=97) ntally
            READ(fp,*,ERR=97) ndata                         
            READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          ENDIF
          IF (debug) WRITE(0,*) '===HEADER COMPLETE==='
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum,binary)
            IF (fluid_ik(icount).NE.0) THEN
              ik = fluid_ik(icount)                          ! Should pull these from .transfer
              ir = fluid_ir(icount)
              IF (imol.EQ.1) THEN                            ! Need check...
                pinmol(ik,ir) = pinmol(ik,ir) + rdum(1)
                pinenm(ik,ir) = pinenm(ik,ir) + rdum(6)
              ELSE
                CALL ER('LoadEireneData_06','IMOL out of bounds, '//
     .                  'unexpected this is...',*99)
              ENDIF
            ENDIF
          ENDDO
          IF (debug) WRITE(0,*) '===DONE==='
c       ----------------------------------------------------------------
        ELSEIF ((binary.AND.code.EQ.-9999932).OR.(.NOT.binary.AND.
     .          buffer(1:22).EQ.'* TEST MOLECULES - SUR')) THEN
          IF (debug) WRITE(0,*) '===TEST MOL: SURFACE FLUXES===',imol
          IF (binary) THEN
            READ(fp  ,ERR=97) ntally
            READ(fp  ,ERR=97) ndata                         
            READ(fp  ,ERR=97) (index(i1),i1=1,ntally)          
          ELSE
            READ(fp,*,ERR=97) ntally
            READ(fp,*,ERR=97) ndata                         
            READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          ENDIF
          icount = 1
          DO WHILE (icount.LE.ndata)
            icount = icount + 1
            CALL NextLine(fp,ntally,iobj,rdum,binary)
            iside = NINT(rdum(1))
            in = wall_in(iside,iobj)
            IF (in.EQ.0) THEN
c             Lots of data comes through that's not associated with 
c             the standard DIVIMP neutral wall, ignore for now...              
c              IF (tri(iobj)%sideindex(4,iside).NE.0) wall_ignored=.TRUE.
              wall_ignored=.TRUE.
              CYCLE
            ENDIF
            IF (imol.EQ.1) THEN                             ! Only for D2...
              IF (tetrahedrons) THEN
                isrf = obj(iobj)%iside(iside)
c               A bit of a hack here, to get the next flux through a surface -- EIRENE stores each surface twice, i.e. one for each 
c               of the associated tetrahedrons, unlike here, where I only store one surface and use a -ve index when linking to
c               it from the "backside", i.e. reverse vertex orientation:              
                fact = SIGN(1.0,REAL(isrf))              
              ELSE
                fact = 1.0
              ENDIF
              tflux(in,7 ) = tflux(in,7 ) + fact * rdum(2)         ! Incident molecule particle flux (s-1)
              tflux(in,8 ) = tflux(in,8 ) + fact * rdum(3)         ! Incident molecule energy flux   (eV s-1)
              tflux(in,13) = tflux(in,13) + fact * SUM(rdum(4:7))  ! Emitted  molecule particle flux (s-1)
              IF (tetrahedrons) THEN
                isrf = ABS(isrf)
                sflux(isrf,1) = 1.0
                sflux(isrf,3) = sflux(isrf,3) +
     .                          fact * (rdum(2) - SUM(rdum(4:7)))

                IF (in.EQ.13) THEN
                  net2 = net2 +  fact * (rdum(2) - SUM(rdum(4:7)))
c                  write(0,'(A,1P,2E12.4,0P,I6,1P,2E12.4,0P)') 
c     .              'net 13 mol:',net2,fact,isrf,rdum(2),SUM(rdum(4:7))
                ENDIF

              ENDIF
              cvesm(in) = 1
c              IF (in.EQ.12) THEN
c                write(0,*) '12: ',iside,iobj,rdum(2),
c     .                            obj(iobj)%iside(iside)
c              ENDIF 
            ENDIF
          ENDDO
          IF (debug) WRITE(0,*) '===DONE==='
c       ----------------------------------------------------------------
        ELSEIF ((binary.AND.code.EQ.-9999941).OR.(.NOT.binary.AND.
     .          buffer(1:17).EQ.'* TEST IONS - VOL')) THEN
          iion = iion + 1
          IF (debug) WRITE(0,*) '===TEST IONS: VOLUME TALLIES===',iion
          IF (binary) THEN
            READ(fp  ,ERR=97) ntally
            READ(fp  ,ERR=97) ndata                         
            READ(fp  ,ERR=97) (index(i1),i1=1,ntally)          
          ELSE
            READ(fp,*,ERR=97) ntally
            READ(fp,*,ERR=97) ndata                         
            READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          ENDIF
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum,binary)
          ENDDO
          IF (debug) WRITE(0,*) '===DONE (NO DATA STORED)==='
c       ----------------------------------------------------------------
        ELSEIF ((binary.AND.code.EQ.-9999942).OR.(.NOT.binary.AND.
     .          buffer(1:17).EQ.'* TEST IONS - SUR')) THEN
          IF (debug) WRITE(0,*) '===TEST IONS: SURFACE FLUXES===',iion
          IF (binary) THEN
            READ(fp  ,ERR=97) ntally
            READ(fp  ,ERR=97) ndata                         
            READ(fp  ,ERR=97) (index(i1),i1=1,ntally)          
          ELSE
            READ(fp,*,ERR=97) ntally
            READ(fp,*,ERR=97) ndata                         
            READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          ENDIF
          icount = 1
          DO WHILE (icount.LE.ndata)
            icount = icount + 1
            CALL NextLine(fp,ntally,iobj,rdum,binary)
            iside = NINT(rdum(1))
            in = wall_in(iside,iobj)
            IF (in.EQ.0) THEN
c             Lots of data comes through that's not associated with 
c             the standard DIVIMP neutral wall, ignore for now...              
c              IF (tri(iobj)%sideindex(4,iside).NE.0) wall_ignored=.TRUE.
              wall_ignored=.TRUE.
              CYCLE
            ENDIF
            IF (iion.EQ.1) THEN                 
            ELSE
              CALL ER('LoadEireneData_06','IION.GT.1, unexpected '//
     .                'this is...',*99)
            ENDIF
          ENDDO
          IF (debug) WRITE(0,*) '===DONE (NO DATA STORED)==='
c       ----------------------------------------------------------------
        ELSEIF ((binary.AND.code.EQ.-9999951).OR.(.NOT.binary.AND.
     .          buffer(1:20).EQ.'* TEST PHOTONS - VOL')) THEN
c       ----------------------------------------------------------------
        ELSEIF ((binary.AND.code.EQ.-9999961).OR.(.NOT.binary.AND.
     .          buffer(1:14).EQ.'* LINE EMISSIO')) THEN
          IF (debug) WRITE(0,*) '===LINE EMISSION==='
          ilin = ilin + 1   
          IF (binary) THEN
            READ(fp  ,ERR=97) ntally
            READ(fp  ,ERR=97) ndata   
            READ(fp  ,ERR=97) (index(i1),i1=1,ntally)          
          ELSE
            READ(fp,*,ERR=97) ntally
            READ(fp,*,ERR=97) ndata   
            READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          ENDIF
          icount = 0

          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum,binary)

            IF (fluid_ik(icount).NE.0) THEN
              ik = fluid_ik(icount)                          ! Should pull these from .transfer
              ir = fluid_ir(icount)
               
              IF     (ilin.EQ.1) THEN
                pinline(ik,ir,1:5,H_BALPHA)=pinline(ik,ir,1:5,H_BALPHA)+ 
     .                                      rdum(1:5)
                pinline(ik,ir,6  ,H_BALPHA)=pinline(ik,ir,6  ,H_BALPHA)+ 
     .                                      rdum(7)
              ELSEIF (ilin.EQ.2) THEN
                pinline(ik,ir,1:6,H_BGAMMA)=pinline(ik,ir,1:6,H_BGAMMA)+ 
     .                                      rdum(1:6)
              ELSEIF (ilin.EQ.3) THEN
                pinline(ik,ir,1:6,H_BBETA )=pinline(ik,ir,1:6,H_BBETA) + 
     .                                      rdum(1:6)
              ENDIF  
            ENDIF
          ENDDO
          IF (debug) WRITE(0,*) '===DONE==='
c       ----------------------------------------------------------------
        ELSEIF ((binary.AND.code.EQ.-9999971).OR.(.NOT.binary.AND.
     .          buffer(1:6 ).EQ.'* MISC')) THEN
c...      Check volumes:
c       ----------------------------------------------------------------
        ELSEIF ((binary.AND.code.EQ.-9999972).OR.(.NOT.binary.AND.
     .          buffer(1:18).EQ.'* PARTICLE SOURCES')) THEN
          IF (debug) WRITE(0,*) '===PARTICLE SORUCES==='
          IF (binary) THEN
            READ(fp  ) i1
          ELSE
            READ(fp,*) i1
          ENDIF

          i2 = 0
          IF (opt_eir%gas_only.EQ.1) i2 = nstrata - gas_n
c          IF (opt_eir%gas_only.EQ.1) i2 = nstrata - 4
          IF ((i1.NE.nstrata  -i2.AND.opt_eir%ntime.EQ.0).OR.
     .        (i1.NE.nstrata+1-i2.AND.opt_eir%ntime.NE.0))
     .      CALL ER('LoadEireneData_06','NSTRATA invalid',*99)

          DO i1 = 1, nstrata
            IF (opt_eir%gas_only.EQ.1.AND.
     .          ((i1.LT.nstrata      .AND.opt_eir%ntime.EQ.0).OR.
     .           (i1.LE.nstrata-gas_n.AND.opt_eir%ntime.NE.0))) THEN
c     .           (i1.LE.nstrata-4.AND.opt_eir%ntime.NE.0))) THEN
              write(0,*) 'skipping strata',i1,nstrata
              strata(i1)%ipanu  = -1
              strata(i1)%fluxt  = -1.0
              strata(i1)%ptrash = -0.01
              strata(i1)%etrash = -0.01
            ELSE
              IF (binary) THEN
                READ(fp  ) i2,strata(i1)%ipanu ,strata(i1)%fluxt ,
     .                        strata(i1)%ptrash,strata(i1)%etrash
              ELSE
                READ(fp,*) i2,strata(i1)%ipanu ,strata(i1)%fluxt ,
     .                        strata(i1)%ptrash,strata(i1)%etrash
              ENDIF
            ENDIF
          ENDDO
c       ----------------------------------------------------------------
        ELSEIF ((binary.AND.code.EQ.-9999973).OR.(.NOT.binary.AND.
     .          buffer(1:13).EQ.'* PUMPED FLUX')) THEN
          IF (debug) WRITE(0,*) '===PUMPED FLUX==='
          DO WHILE (.TRUE.) 
            IF (binary) THEN
              READ(fp      ) isur,species(1:6),amps
                write(0,*) '>>>>>>>>>'//species(1:6)//'<',amps
              IF (isur.EQ.-9999999) THEN
                goodeof = .TRUE.
                EXIT
              ENDIF
            ELSE
              READ(fp,'(A256)',END=97,ERR=97) buffer
              IF (buffer(1:6).EQ.'* DONE') THEN
                goodeof = .TRUE.
                EXIT
              ENDIF
              READ(buffer,*) isur,species,amps
            ENDIF
            IF     (species(1:6).EQ.'D     '.OR.                   ! Need to do something to record the 
     .              species(1:6).EQ.'D(N=1)') THEN                 ! flux pumped by the core boundary
              pflux = amps / ECH 
              hescpd  = hescpd + pflux
              IF (isur.NE.-core_boundary) hescal = hescal + pflux    
            ELSEIF (species(1:6).EQ.'D2    ') THEN
              pflux = amps / ECH * 2.0
              hescpd  = hescpd + pflux
              IF (isur.NE.-core_boundary) hescal = hescal + pflux
            ELSE
              CALL WN('LoadEireneData','Unknown particle type')
              WRITE(0,*) '  LABEL= '//TRIM(species)
            ENDIF          
          ENDDO
          IF (debug) WRITE(0,*) '===DONE==='
c       ----------------------------------------------------------------
        ELSEIF ((binary.AND.code.EQ.-9999974).OR.(.NOT.binary.AND.
     .          buffer(1:16).EQ.'* ITERATION DATA')) THEN
          IF (debug) WRITE(0,*) '===ITERATION DATA==='
          IF (binary) THEN
            READ(fp  ,END=97,ERR=97) i1
          ELSE
            READ(fp,*,END=97,ERR=97) i1
          ENDIF
          WRITE(0,*) 'i1 nhistory',i1
          nhistory = nhistory + 1
          DO i2 = 1, i1
            IF (binary) THEN
              READ(fp  ,END=97,ERR=97) i3
            ELSE
              READ(fp,*,END=97,ERR=97)
              READ(fp,*,END=97,ERR=97)
              READ(fp,*,END=97,ERR=97) i3
            ENDIF
            history(nhistory)%iiter   = iitersol
            history(nhistory)%ngauge  = i1
            history(nhistory)%nstrata = i3
c            WRITE(0,*) 'i2,3',i2,i3
            IF (binary) THEN
            ELSE
              READ(fp,*,END=97,ERR=97) 
            ENDIF
c            READ(fp,*,END=97,ERR=97) 
            DO i4 = 1, i3
              IF (binary) THEN
                READ(fp  ,END=97,ERR=97) 
     .            idum,
     .            history(nhistory)%gauge_vol       (   i2),   ! [m-3]
     .            history(nhistory)%gauge_p_atm     (i4,i2),   ! [mTorr]
     .            history(nhistory)%gauge_parden_atm(i4,i2),   ! [particles m-3]
     .            history(nhistory)%gauge_egyden_atm(i4,i2),   ! [eV m-3]
     .            history(nhistory)%gauge_p_mol     (i4,i2), 
     .            history(nhistory)%gauge_parden_mol(i4,i2), 
     .            history(nhistory)%gauge_egyden_mol(i4,i2) 
              ELSE
                READ(fp,*,END=97,ERR=97) 
     .            idum,
     .            history(nhistory)%gauge_vol       (   i2),   ! [m-3]
     .            history(nhistory)%gauge_p_atm     (i4,i2),   ! [mTorr]
     .            history(nhistory)%gauge_parden_atm(i4,i2),   ! [particles m-3]
     .            history(nhistory)%gauge_egyden_atm(i4,i2),   ! [eV m-3]
     .            history(nhistory)%gauge_p_mol     (i4,i2), 
     .            history(nhistory)%gauge_parden_mol(i4,i2), 
     .            history(nhistory)%gauge_egyden_mol(i4,i2) 
              ENDIF
            ENDDO
          ENDDO
c       ----------------------------------------------------------------
        ELSEIF (.NOT.binary.AND.buffer(1:1 ).EQ.'*') THEN
c       ----------------------------------------------------------------
        ELSE
        ENDIF

      ENDDO
 10   CONTINUE

      CLOSE (fp)

c...  Relaxation into OSM arrays only?

c...  Normalize volume quantities (need to be careful vis-a-vis relaxation):
      sumion = 0.0
      zioniz = 0.0
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        DO ik = 1, nks(ir)
          norm = 1.0 / kvols(ik,ir)  ! KVOLS scaled by EIRTORFRAC in tau.f

c...      Volume normalization:
          tdata(ik,ir,1:5) = tdata(ik,ir,1:5) * norm
c...      Linear relaxation:
          pinion(ik,ir) = (1.0-frac)*pinion(ik,ir) + frac*tdata(ik,ir,1)
          pinrec(ik,ir) = (1.0-frac)*pinrec(ik,ir) + frac*tdata(ik,ir,2)
          pinmp (ik,ir) = (1.0-frac)*pinmp (ik,ir) + frac*tdata(ik,ir,3)
          pinqi (ik,ir) = (1.0-frac)*pinqi (ik,ir) + frac*tdata(ik,ir,4)
          pinqe (ik,ir) = (1.0-frac)*pinqe (ik,ir) + frac*tdata(ik,ir,5)
c...      Not relaxed, volume normalize:
          pinatom(ik,ir) = pinatom(ik,ir) * norm  
          pinena (ik,ir) = pinena (ik,ir) * norm  
          pinmol (ik,ir) = pinmol (ik,ir) * norm  
          pinenm (ik,ir) = pinenm (ik,ir) * norm  
          pinline(ik,ir,1:6,H_BALPHA) = pinline(ik,ir,1:6,H_BALPHA) * !/
     .                                  norm  
          pinline(ik,ir,1:6,H_BGAMMA) = pinline(ik,ir,1:6,H_BGAMMA) * !/
     .                                  norm  
          pinline(ik,ir,1:6,H_BBETA)  = pinline(ik,ir,1:6,H_BBETA)  * !/
     .                                  norm  
c...
          pinalpha(ik,ir) = pinline(ik,ir,6,H_BALPHA)
c...  
          IF (ilspt.GT.0) THEN
            eirdat(ik,ir,:) = eirdat(ik,ir,:) * norm
            pinenz(ik,ir  ) = pinenz(ik,ir  ) * norm

            pinionz(ik,ir) = eirdat(ik,ir,1)  ! Impurity ionisation rate
            pinz0  (ik,ir) = eirdat(ik,ir,2)  ! Density of neutral impurity particles

            zioniz = zioniz + pinionz(ik,ir) * karea2(ik,ir)               
          ENDIF
c...      
          IF (pinatom(ik,ir).GT.1.0E-10) THEN
            pinena(ik,ir) = pinena(ik,ir) / pinatom(ik,ir)
          ELSE
            pinena(ik,ir) = 0.0
          ENDIF
          IF (pinmol (ik,ir).GT.1.0E-10) THEN
            pinenm(ik,ir) = pinenm(ik,ir) / pinmol (ik,ir)
          ELSE
            pinenm(ik,ir) = 0.0
          ENDIF     
          IF (pinz0  (ik,ir).GT.1.0E-10) THEN
            pinenz(ik,ir) = pinenz(ik,ir) / pinz0  (ik,ir)
          ELSE
            pinenz(ik,ir) = 0.0
          ENDIF

          sumion = sumion + pinion(ik,ir) * kvols(ik,ir) * eirtorfrac
        ENDDO
      ENDDO

      WRITE(PINOUT,*) 'SUMION:',sumion

c...  Assign wall fluxes (based on code in ReadWallFlux in 
c     setup.f, for EIRENE99):

c     FLUXHW - FLUX OF HYDROGEN (ATOMS AND MOLECULES) TO THE WALL
c     FLXHW2 - FLUX OF HYDROGEN (ATOMS AND IONS) TO THE WALL
c     FLXHW3 - FLUX OF IMPURITIES SPUTTERED FROM THE WALL (N/A)
c     FLXHW4 - FLUX OF IMPURITIES REDEPOSITED ONTO THE WALL (N/A)  --- *HACK* AVERAGE IMPURITY LAUNCH ENERGY
c     FLXHW5 - AVERAGE ENERGY OF ATOMS HITTING THE WALL (EV)
c     FLXHW6 - FLUX OF HYDROGEN ATOMS TO THE WALL
c     FLXHW7 - AVERAGE ENERGY OF MOLECULES HITTING THE WALL (eV)
c     FLXHW8 - EIRENE REPORTED HYDROGEN ION FLUXES TO THE WALL 
c
      fluxhw = 0.0
      flxhw2 = 0.0
      flxhw3 = 0.0
      flxhw4 = 0.0
      flxhw5 = 0.0
      flxhw6 = 0.0
      flxhw7 = 0.0
      flxhw8 = 0.0

      DO iw = 1, nvesm+nvesp
        cir = 2.0 * PI * 0.5 * (rvesm(iw,1) + rvesm(iw,2))
        len = SQRT((rvesm(iw,1) - rvesm(iw,2))**2.0 +
     .             (zvesm(iw,1) - zvesm(iw,2))**2.0)
        area = cir * len
c        area = 1.0 / ECH
        fluxhw(iw) = (tflux(iw,1) + tflux(iw,7)) / area    ! *** BUG? *** since mixing D and D2?
        flxhw2(iw) = (tflux(iw,1) + tflux(iw,3)) / area  
        flxhw3(iw) = tflux(iw,5) / area
        IF (tflux(iw,5).NE.0.0) flxhw4(iw) = tflux(iw,6) / tflux(iw,5)  ! *** HACK AVERAGE IMPURITY INJECTION ENERGY ***
        IF (tflux(iw,1).NE.0.0) flxhw5(iw) = tflux(iw,2) / tflux(iw,1)  ! flxhw5(iw) = tflux(iw,2) / (tflux(iw,1) + 1.0E-10)
        flxhw6(iw) = tflux(iw,1) / area 
        IF (tflux(iw,7).NE.0.0) flxhw7(iw) = tflux(iw,8) / tflux(iw,7)  
c...    If this is fixed on the EIRENE side so that it is no longer
c       statistical, then the use of FLXHW8 in the CALC_TARGFLUXDATA 
c       routine can be removed:
        flxhw8(iw) = tflux(iw,3) / area  !  * 1.602E-19  For hand checks...
      ENDDO
c...  Check that flux data was assigned to every wall/target segment:
      warning_message = .FALSE.
      DO iw = 1, nvesm+nvesp
        IF (cvesm(iw).EQ.0) THEN
          warning_message = .TRUE.
          WRITE(eirfp,'(A,3I4)') 'WARNING: NO WALL FLUX DATA ASS'//
     .                           'IGNED FOR SEGMENT ',iw,nvesm,nvesp
        ENDIF
      ENDDO
      IF (warning_message) THEN
        WRITE(0,*) 
        WRITE(0,*) '*********************************************'
        WRITE(0,*) '* NO FLUX DATA FOUND FOR SOME WALL SEGMENTS *'
        WRITE(0,*) '*********************************************'
        WRITE(0,*) 
      ENDIF
c...  Store data in allocatable structures:
      IF (.NOT.ALLOCATED(wall_flx)) THEN
        wall_n       = nvesm + nvesp
        wall_nlaunch = 0
        ALLOCATE(wall_flx(wall_n))
        call init_wall_flx(wall_n)
        ! jdemod - created initialization routine in module since
        !          this can also be allocated in neut.f for some reason
        !
        !DO i = 1, wall_n
        !  wall_flx(i)%in_par_blk = 0.0
        !  wall_flx(i)%in_par_atm = 0.0
        !  wall_flx(i)%in_par_mol = 0.0
        !  wall_flx(i)%em_par_mol = 0.0
        !  wall_flx(i)%in_ene_blk = 0.0
        !  wall_flx(i)%in_ene_atm = 0.0
        !  wall_flx(i)%in_ene_mol = 0.0
        !  wall_flx(i)%em_par_atm = 0.0
        !  wall_flx(i)%em_ene_atm = 0.0
        !  wall_flx(i)%launch     = 0.0        
        !  wall_flx(i)%prompt     = 0.0        
        !ENDDO
      ENDIF
      DO i = 1, wall_n
        cir = 2.0 * PI * 0.5 * (rvesm(i,1) + rvesm(i,2))
        len = SQRT((rvesm(i,1) - rvesm(i,2))**2.0 +
     .             (zvesm(i,1) - zvesm(i,2))**2.0)
        area = cir * len
c        area = 1.0 / ECH
        wall_flx(i)%length = len
        wall_flx(i)%area   = area

c       The first index is for the species index for each species type, i.e
c       atom (_atm) type number X, where X for atoms is usually equal to 1 (D) or 
c       2 (D and a sputtered impurity species).
c
c       The second index is for the species type of the particle source that's
c       responsible for producing the incident/emitted particle, following
c       the EIRENE convention:
c
c          0-total,1-bulk particles,2-test atoms,3-test mol.,4-test ions,5-photons

        wall_flx(i)%in_par_blk(1,0) = flxhw8(i)
        wall_flx(i)%in_par_atm(1,0) = flxhw6(i)
        wall_flx(i)%em_par_atm(1,0) = tflux(i,14) / area
        wall_flx(i)%in_par_mol(1,0) = fluxhw(i) - flxhw6(i)
        wall_flx(i)%em_par_mol(1,0) = tflux(i,13) / area

c        wall_flx(i)%in_par_mol(1,0) = fluxhw(i)  ! *** BUG ***  FLUXHW is not the correct quantitiy 10/03/2011 -SL

        IF (tflux(i,3).NE.0.0) 
     .    wall_flx(i)%in_ene_blk(1,0) = tflux(i,4) / tflux(i,3)

        wall_flx(i)%in_ene_atm(1,0) = flxhw5(i)
        wall_flx(i)%in_ene_mol(1,0) = flxhw7(i)

        wall_flx(i)%in_par_atm(2,0) = tflux(i,11) / area         ! Total impurity atom flux to the surface
        IF (tflux(i,11).NE.0.0) 
     .  wall_flx(i)%in_ene_atm(2,0) = tflux(i,12) / tflux(i,11)  ! Average energy of the impurity atoms hitting the surface

        wall_flx(i)%em_par_atm(2,1) = tflux(i,9) / area                  ! Impurity atom influx from bulk ions
        wall_flx(i)%em_par_atm(2,2) = tflux(i,5) / area                  ! Impurity atom influx from test atoms
        IF (tflux(i,9).NE.0.0) wall_flx(i)%em_ene_atm(2,1)=tflux(i,10) /
     .                                                     tflux(i,9 )
        IF (tflux(i,5).NE.0.0) wall_flx(i)%em_ene_atm(2,2)=tflux(i,6 ) / 
     .                                                     tflux(i,5 )
      ENDDO
c...  Not sure what PINCOR is about at the moment, so issue a warning
c     if it is not unity:
      IF (pincor.NE.1.0) THEN
        CALL WN('LoadEireneData_06','PINCOR not unity, setting to 1.0')      
        pincor = 1.0
      ENDIF
c...  Not sure the reason at the moment, but this is done in READPIN 
c     and the _PIN arrays are referenced in TAU:
      DO iw = 1, nvesm+nvesp
        fluxhw(iw) = fluxhw(iw) * pincor
        flxhw2(iw) = flxhw2(iw) * pincor
        flxhw3(iw) = flxhw3(iw) * pincor
        flxhw4(iw) = flxhw4(iw) * pincor
        flxhw6(iw) = flxhw6(iw) * pincor
c...    Copy wall flux data:
        fluxhw_pin(iw) = fluxhw(iw)
        flxhw2_pin(iw) = flxhw2(iw)
        flxhw3_pin(iw) = flxhw3(iw)
        flxhw4_pin(iw) = flxhw4(iw)
        flxhw5_pin(iw) = flxhw5(iw)
        flxhw6_pin(iw) = flxhw6(iw)
c...    Copy vessel definitions:
        jvesm_pin(iw)  = jvesm(iw)
        DO in = 1,2
          rvesm_pin(iw,in)  = rvesm(iw,in)
          zvesm_pin(iw,in)  = zvesm(iw,in)
        ENDDO
      ENDDO
c...  From READPIN:
      IF ((wlpabs.EQ.2.OR.wlpabs.EQ.3).AND.cgeoopt.NE.-1) 
     .  CALL ER('LoadEireneData_06','Probability data not passed '//
     .                              'from EIRENE',*99)
c...  Also from READPIN?
      hcorr = 1.0  ! volume correction factor
      hval  = 1.0  ! volume correction factor

      CALL OutputEIRENE(67,'WORKING WITH EIRENE06')

c...  Dump EIRENE calculated impurity distribution data:
      IF (ALLOCATED(eirdat)) THEN
        CALL inOpenInterface('idl.eirene_imp',ITF_WRITE)
        CALL inPutData(0.0       ,'IMP_INITIAL_IZ'    ,'N/A')
        CALL inPutData(0.0       ,'IMP_MAX_IZ'        ,'N/A')
        CALL inPutData(REAL(cion),'IMP_Z'             ,'N/A')
        CALL inPutData(crmi      ,'IMP_A'             ,'N/A')
        CALL inPutData(irsep -1  ,'GRID_ISEP'         ,'N/A')  ! TUBE is set to the OSM fluid grid system, where                   
        CALL inPutData(irwall-1  ,'GRID_IPFZ'         ,'N/A')  ! the boundary rings are not present
        CALL inPutData(eirtorfrac,'TOROIDAL_FRACTION' ,'N/A')  
        DO ir = 2, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          ike = nks(ir)
          IF (ir.LT.irsep) ike = ike - 1
          tube = ir - 1                      
          IF (ir.GT.irwall) tube = tube - 2  
          DO ik = 1, ike
            CALL inPutData(-1  ,'INDEX','N/A')                     
            CALL inPutData(ik  ,'POS'  ,'N/A')                     
            CALL inPutData(tube,'TUBE' ,'N/A')  
            CALL inPutData(kss(ik,ir),'S','m')
            CALL inPutData(kps(ik,ir),'P','m')
            CALL inPutData(kvols(ik,ir),'VOLUME','m-3')
            CALL inPutData(eirdat(ik,ir,2),'IMP_DENS_00' ,'m-3 s-1')  
            CALL inPutData(eirdat(ik,ir,1),'IMP_IONIZ_00','m-3 s-1')  
          ENDDO
        ENDDO
        CALL inCloseInterface
        DEALLOCATE(eirdat)
      ENDIF

c...  Dump EIRENE iteration data:
      CALL inOpenInterface('idl.eirene_history',ITF_WRITE)
      CALL inPutData(opt_eir%time0*1.0E-6,'TIME0','s')
      CALL inPutData(time0        *1.0E-6,'TIME' ,'s')
      DO i1 = 1, nhistory
       DO i2 = 1, history(i1)%ngauge
        CALL inPutData(history(i1)%gauge_vol(i2),'VOLUME','m-3')
        DO i3 = 1, history(i1)%nstrata
         CALL inPutData(opt_eir%gauge_ind(1,i2),'GAUGE_I1' ,'N/A')        !
         CALL inPutData(opt_eir%gauge_ind(2,i2),'GAUGE_I2' ,'N/A')        !
         CALL inPutData(i3,'STRATA' ,'N/A')                !
         CALL inPutData(i1,'HISTORY','N/A')                !
         CALL inPutData(i2,'GAUGE'  ,'N/A')                !
         CALL inPutData(opt_eir%gauge_pos(1,i2),'GAUGE_X'  ,'N/A')        !
         CALL inPutData(opt_eir%gauge_pos(2,i2),'GAUGE_Y'  ,'N/A')        !
         CALL inPutData(opt_eir%gauge_pos(3,i2),'GAUGE_Z'  ,'N/A')        !
         CALL inPutData(opt_eir%gauge_pos(4,i2),'GAUGE_PHI','N/A')        !
         CALL inPutData(history(i1)%iiter,'FLUID_ITERATION','N/A')        !
         rdum(1) = history(i1)%gauge_p_atm(i3,i2) / 7.502  ! from 101.3 Pa = 760 mTorr
         rdum(2) = history(i1)%gauge_p_mol(i3,i2) / 7.502  
         rdum(3) = history(i1)%gauge_egyden_atm(i3,i2) /
     .            (history(i1)%gauge_parden_atm(i3,i2) + 1.0E-10)
         rdum(4) = history(i1)%gauge_egyden_mol(i3,i2) /
     .            (history(i1)%gauge_parden_mol(i3,i2) + 1.0E-10)
         rdum(5) = history(i1)%gauge_parden_atm(i3,i2)
         rdum(6) = history(i1)%gauge_parden_mol(i3,i2)
         CALL inPutData(history(i1)%gauge_p_atm(i3,i2),'P1_ATM','mTorr')
         CALL inPutData(history(i1)%gauge_p_mol(i3,i2),'P1_MOL','mTorr')
         CALL inPutData(rdum(1)                       ,'P2_ATM','Pa')
         CALL inPutData(rdum(2)                       ,'P2_MOL','Pa')
         CALL inPutData(rdum(3)                       ,'EAVG_ATM','eV')
         CALL inPutData(rdum(4)                       ,'EAVG_MOL','eV')
         CALL inPutData(rdum(5)                       ,'DENS_ATM','m-3')
         CALL inPutData(rdum(6)                       ,'DENS_MOL','m-3')
        ENDDO
       ENDDO
      ENDDO
      CALL inCloseInterface

c...  Saving wall flux data:
      IF (ALLOCATED(wall_flx)) THEN
        CALL inOpenInterface('idl.eirene_flux_wall',ITF_WRITE)
        CALL inPutData(wall_n      ,'N_SEGMENTS','N/A')
        CALL inPutData(wall_nlaunch,'N_LAUNCH'  ,'N/A')
        CALL inPutData(MAXNLAUNCH  ,'MAXNLAUNCH','N/A')
        CALL inPutData(MAXNBLK     ,'MAXNBLK'   ,'N/A')
        CALL inPutData(MAXNATM     ,'MAXNATM'   ,'N/A')
        CALL inPutData(MAXNMOL     ,'MAXNMOL'   ,'N/A')
        CALL inPutData(MAXNION     ,'MAXNION'   ,'N/A')
        CALL inPutData(MAXNPHO     ,'MAXNPHO'   ,'N/A')
        CALL inPutData(MAXNSRC     ,'MAXNSRC'   ,'N/A')
        CALL inPutData(SUM(wall_flx(:)%em_par_atm(2,1) *
     .                     wall_flx(:)%length),'TOT_EM_IMP_1','s-1')
        CALL inPutData(SUM(wall_flx(:)%em_par_atm(2,2) *
     .                     wall_flx(:)%length),'TOT_EM_IMP_2','s-1')
        CALL inPutData(SUM(wall_flx(:)%em_par_atm(2,1) * ECH *
     .                     wall_flx(:)%area  ),'TOT_EM_IMP_1','s-1')
        CALL inPutData(SUM(wall_flx(:)%em_par_atm(2,2) * ECH * 
     .                     wall_flx(:)%area  ),'TOT_EM_IMP_2','s-1')
        CALL inPutData(wall_flx(:)%length,'LENGTH','m' )
        CALL inPutData(wall_flx(:)%area  ,'AREA  ','m2')
        DO i = 1, MAXNBLK
          WRITE(tag,'(A,I1,A,10X)') 'IN_PAR_BLK_',i,'_0'
          CALL inPutData(wall_flx(:)%in_par_blk(i,0),tag,'m-2 s-1')
          WRITE(tag,'(A,I1,A,10X)') 'IN_ENE_BLK_',i,'_0'
          CALL inPutData(wall_flx(:)%in_ene_blk(i,0),tag,'eV')
        ENDDO
        DO i = 1, MAXNATM
          WRITE(tag,'(A,I1,A,10X)') 'IN_ENE_ATM_',i,'_0'
          CALL inPutData(wall_flx(:)%in_ene_atm(i,0),tag,'eV')
          WRITE(tag,'(A,I1,A,10X)') 'IN_PAR_ATM_',i,'_0'
          CALL inPutData(wall_flx(:)%in_par_atm(i,0),tag,'m-2 s-1')
          WRITE(tag,'(A,I1,A,10X)') 'EM_PAR_ATM_',i,'_0'
          CALL inPutData(wall_flx(:)%em_par_atm(i,0),tag,'m-2 s-1')
        ENDDO
        DO i = 1, MAXNMOL
          WRITE(tag,'(A,I1,A,10X)') 'IN_ENE_MOL_',i,'_0'
          CALL inPutData(wall_flx(:)%in_ene_mol(i,0),tag,'eV')
          WRITE(tag,'(A,I1,A,10X)') 'IN_PAR_MOL_',i,'_0'
          CALL inPutData(wall_flx(:)%in_par_mol(i,0),tag,'m-2 s-1')
          WRITE(tag,'(A,I1,A,10X)') 'EM_PAR_MOL_',i,'_0'
          CALL inPutData(wall_flx(:)%em_par_mol(i,0),tag,'m-2 s-1')
          WRITE(tag,'(A,I1,A,10X)') 'NT_PAR_',i,'_0'
          CALL inPutData(2.0*wall_flx(:)%in_par_mol(i,0) + 
     .                       wall_flx(:)%in_par_atm(i,0) -
     .                   2.0*wall_flx(:)%em_par_mol(i,0) -
     .                       wall_flx(:)%em_par_atm(i,0), 
     .                   tag,'m-2 s-1')

          write(0,*) 'total through=',
     .      SUM(2.0*wall_flx(:)%in_par_mol(i,0) + 
     .              wall_flx(:)%in_par_atm(i,0) -
     .          2.0*wall_flx(:)%em_par_mol(i,0) -
     .              wall_flx(:)%em_par_atm(i,0))

        ENDDO
        DO i = 2, 2
          WRITE(tag,'(A,I1,A,10X)') 'EM_PAR_ATM_',i,'_1'
          CALL inPutData(wall_flx(:)%em_par_atm(i,1),tag,'m-2 s-1')
          WRITE(tag,'(A,I1,A,10X)') 'EM_ENE_ATM_',i,'_1'
          CALL inPutData(wall_flx(:)%em_ene_atm(i,1),tag,'eV')
          WRITE(tag,'(A,I1,A,10X)') 'EM_PAR_ATM_',i,'_2'
          CALL inPutData(wall_flx(:)%em_par_atm(i,2),tag,'m-2 s-1')
          WRITE(tag,'(A,I1,A,10X)') 'EM_ENE_ATM_',i,'_2'
          CALL inPutData(wall_flx(:)%em_ene_atm(i,2),tag,'eV')
        ENDDO
        DO i = 1, wall_nlaunch
          WRITE(tag,'(A,I0.2,10X)') 'LAUNCH_',i
          CALL inPutData(wall_flx(:)%launch(i),tag,'???')
        ENDDO
        CALL inPutData(wall_flx(:)%prompt,'PROMPT_DEP','s-1 m-2')
        CALL inCloseInterface
      ENDIF


c...  Dump EIRENE surface flux data:
      IF (tetrahedrons) THEN 

        ALLOCATE(vmap(2,nvtx))
        vmap = 0

        CALL inOpenInterface('idl.tet_flux',ITF_WRITE)

        CALL inPutData(1.1,'version','N/A')

        j = 0
        DO isrf = 1, nsrf
          IF (sflux(isrf,1).EQ.0.0) CYCLE   

c         Calculate the area of the 3D triangle, using the half cross product, from here:
c           http://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
          v1(1) = vtx(1,srf(isrf)%ivtx(1)) - vtx(1,srf(isrf)%ivtx(2))
          v1(2) = vtx(2,srf(isrf)%ivtx(1)) - vtx(2,srf(isrf)%ivtx(2))
          v1(3) = vtx(3,srf(isrf)%ivtx(1)) - vtx(3,srf(isrf)%ivtx(2))
          v2(1) = vtx(1,srf(isrf)%ivtx(3)) - vtx(1,srf(isrf)%ivtx(2))
          v2(2) = vtx(2,srf(isrf)%ivtx(3)) - vtx(2,srf(isrf)%ivtx(2))
          v2(3) = vtx(3,srf(isrf)%ivtx(3)) - vtx(3,srf(isrf)%ivtx(2))

          area = SQRT( (v1(2) * v2(3) - v1(3) * v2(2))**2 +
     .                 (v1(3) * v2(1) - v1(1) * v2(3))**2 +
     .                 (v1(1) * v2(2) - v1(2) * v2(1))**2 ) * 0.5

c         Identify the EIRENE surface associated with this triangle:
          iliin = -1
          DO i = 1, nsurface
            IF (surface(i)%type.NE.NON_DEFAULT_STANDARD) CYCLE
            IF (surface(i)%num.EQ.srf(isrf)%index(IND_SURFACE)) THEN
              iliin = surface(i)%iliin 
            ENDIF
          ENDDO          
          CALL inPutData(isrf,'ISRF','N/A')
          DO i = 1, 3
            ivtx = srf(isrf)%ivtx(i)
            IF (vmap(1,ivtx).EQ.0) THEN
              j = j + 1
              vmap(1,ivtx) = j
              vmap(2,j   ) = ivtx
            ENDIF
            WRITE(tag,'(A,I0.1)') 'V',i
            CALL inPutData(vmap(1,ivtx),TRIM(tag),'N/A')                  
          ENDDO
          CALL inPutData(srf(isrf)%index(IND_WALL_STD),'IN_STD','N/A')
          CALL inPutData(srf(isrf)%index(IND_WALL_ADD),'IN_ADD','N/A')
          CALL inPutData(iliin,'ILIIN','N/A')
          CALL inPutData(area ,'AREA' ,'m-2')
          CALL inPutData(sflux(isrf,2)/area,'NT_PAR_ATM_1_0','m-2 s-1')
          CALL inPutData(sflux(isrf,3)/area,'NT_PAR_MOL_1_0','m-2 s-1')
          CALL inPutData(sflux(isrf,4)/area,'NT_PAR_IMP_1_0','m-2 s-1')
        ENDDO

c       Write the referenced vertices for this group of surfaces:
        DO i = 1, j
          CALL inPutData(vtx(1,vmap(2,i)),'X','m')
          CALL inPutData(vtx(2,vmap(2,i)),'Y','m')
          CALL inPutData(vtx(3,vmap(2,i)),'Z','m')         
        ENDDO


        CALL inCloseInterface

        DEALLOCATE(vmap)
      ENDIF

c...  Need to save data again since Dalpha has been loaded into OBJ:
c      CALL SaveGeometryData('tetrahedrons.raw')

      write(0,*) 'clearing objects',nobj,nsrf,nvtx

c...  Clear arrays: 
      DEALLOCATE(fluid_ik)
      DEALLOCATE(fluid_ir)
      DEALLOCATE(wall_in)
      DEALLOCATE(tdata)
      DEALLOCATE(tflux)
      IF (ALLOCATED(sflux)) DEALLOCATE(sflux)
      IF (ALLOCATED(tri  )) DEALLOCATE(tri  ) ! Move earlier and pass index mapping to/from Eirene
      IF (ALLOCATED(obj  )) DEALLOCATE(obj  ) ! *same* but need to make more selective
      IF (ALLOCATED(srf  )) DEALLOCATE(srf  ) ! at some point since even fluid objects will be stored here... 
      IF (ALLOCATED(vtx  )) DEALLOCATE(vtx  )

c...  Output results of confidence checks:
      IF (.NOT.goodeof)  CALL ER('LoadEireneData','Problem with '//
     .                           'eirene.transfer file',*99)
      IF (wall_ignored) THEN
        WRITE(0,*) 
        WRITE(0,*) '----------------------------------------------'
        WRITE(0,*) '  FLUXES TO NON-STANDARD DIVIMP WALL IGNORED  '
        WRITE(0,*) '----------------------------------------------'
        WRITE(0,*) 
      ENDIF

      RETURN
 97   CALL ER('LoadEireneData_06','Problem reading data file',*99)
 98   CALL ER('LoadEireneData_06','Data file not found',*99)
 99   WRITE(0,*) 'IOBJ =',icount
      WRITE(0,*) 'ISIDE=',iside
      WRITE(0,*) '  FILENAME= "'//TRIM(fname)//'"'
      WRITE(0,*) 'NSTRATA=',nstrata
      WRITE(0,*) 'I1     =',i1
      STOP
      END
c
c ======================================================================
c


