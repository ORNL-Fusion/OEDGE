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
      USE mod_sol28_global
      USE mod_geometry
      USE mod_filament
      IMPLICIT none

      include 'params'
      include 'cgeom'
      include 'comtor'
      INCLUDE 'slcom'

      INTEGER, INTENT(IN) :: iitersol

      INTEGER ik,ir,in1,in2,i1,id,ik1,status
      LOGICAL saved_triangles,output

      REAL*8 t  ! *** TEMP ***

      DATA t /0.0D0/
      DATA saved_triangles /.FALSE./
      SAVE

      IF (citersol.GT.0) THEN
        WRITE(0,*) 
        WRITE(0,*) '---------------------------------'
        WRITE(0,*) iitersol,SNGL(t),nfilament
        WRITE(0,*) '---------------------------------'
        WRITE(0,*) 
        IF (t.EQ.0.0D0) CALL DefineFilaments 
        CALL SetupFilaments(t)
        t = t + 10.0D-06
      ELSE
        IF (t.EQ.0.0D0) CALL DefineFilaments 
        CALL SetupFilaments(t)
      ENDIF

c      CALL DefineFilaments
c      t = 1.0D0
c      DO i1 = 1, 10
c        CALL SetupFilaments(t)      
c        t = t + 10.0D-6
c      ENDDO
c      STOP 'asshole'

      helium = .TRUE.

      output = .FALSE.
      opt%pin_data = .TRUE.
      tetrahedrons = .FALSE.

      SELECTCASE (eirgeom)
        CASE (-1)
          CALL ER('WriteEireneFiles_06','Bad EIRGEOM',*99)
        CASE ( 3)
          tetrahedrons = .TRUE.
      ENDSELECT        

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
      dtimv = eirdtimv                          ! Time-dependent runs

      ttemp = ctargt                            ! Temperature (K) of target surfaces
      wtemp = cwallt                            ! Temperature (K) of wall surfaces

c...  Surface material, for use with the local reflection models
c     specified in block 3 in the EIRENE input file:
      IF (eirmat1.EQ.1) tmater = 9642.0         ! Target material
      IF (eirmat1.EQ.2) tmater = 1206.0
      IF (eirmat1.EQ.3) tmater = 18474.0
      IF (eirmat1.EQ.4) tmater = 904.0  
      IF (eirmat2.EQ.1) wmater = 9642.0         ! Wall
      IF (eirmat2.EQ.2) wmater = 1206.0
      IF (eirmat2.EQ.3) wmater = 18474.0
      IF (eirmat2.EQ.4) wmater = 904.0

      opacity   = eiropacity
      photons   = eirphoton
      trim_data = eirtrim
      bgk       = eirbgk
      ntorseg   = eirntorseg
      torfrac   = eirtorfrac
      alloc     = eiralloc
      beam      = 0

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

        IF (saved_triangles.AND..NOT.tetrahedrons) THEN
          IF (tetrahedrons) THEN          
            STOP 'need to save plasma / target data with .raw'
          ELSE
            CALL LoadTriangles_06
          ENDIF
          CALL ProcessFluidGrid_06
          CALL AssignPlasmaQuantities_06
          CALL SetupEireneStrata
        ELSE
          IF (output) WRITE(0,*) 'building triangles'
          CALL ALLOC_VERTEX  (10000)   ! NEED TO ADD ACTIVE BOUNDS CHECKING FOR ALL THESE!
          CALL ALLOC_SURFACE (300)
          CALL ALLOC_TRIANGLE(25000)

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
          CALL WritePolyFile_06(eirntri,MAXNRS,eirtri)

c...      Re-builds connection map, taking into account the new triangles:
          CALL ProcessTriangles_06(0)

c...      Dumps trianlges to a binary file, for use with LoadTriangles:
          CALL SaveTriangles_06
          saved_triangles = .TRUE.

c...      Tetrahedrons:
          IF (tetrahedrons) THEN
            CALL ProcessTetrahedrons_06
            CALL SaveGeometryData('tetrahedrons.raw')
c            CALL DumpGrid('BUILDING TETRAHEDRONS')          
          ENDIF

          WRITE(0,*) 'DONE'
        ENDIF


        IF (.TRUE.) CALL LocalGridRefinement


c...    Writes the .points, .sides, .map and .plasma files to be loaded
c       by EIRENE:
        IF (tetrahedrons) THEN
          CALL WriteEireneObjects
        ELSE        
          CALL WriteEireneTriangles
        ENDIF

        IF (eirdata.EQ.1) CALL WriteEireneInputFile_06
      ENDIF

c...  Fluid grid data no longer required:
      CALL DEALLOC_CELL

c...  Need this to generate celldata.dat for triangle grid, which
c     is required for finding the B-field everywhere in the vessel:
c      CALL ProcessTriangles

c      IF (tetrahedrons) CALL DumpGrid('BUILDING TETRAHEDRONS')  

c...  TMP: can be removed...
c     (EDGE2D: not required, OSM specific)
c      IF     (s28cfpdrft.GE.1) THEN
c        WRITE(0,*) 'CALCULATING RADIAL FLUX AFTER RETURNING FROM EIRENE'
c        CALL CalcRadialDrift(-1)
c        CALL CalcRadialDrift(-2)
c      ELSEIF (s28cfpdrft.EQ.-1) THEN
cc...    Turning it off:
c        WRITE(0,*) 'TURING OFF RADIAL FLUX AFTER RETURNING FROM EIRENE'
c        DO ir = 1, nrs
c          DO ik = 1, nks(ir)
c            osmcfpflx(ik,ir,2) = (1.0 - rel_frac) * osmcfpflx(ik,ir,2) 
c          ENDDO
c        ENDDO
c        CALL CalcRadialDrift(-2)
c      ENDIF

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: DefineEireneSurfaces
c
      SUBROUTINE DefineEireneSurfaces_06
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none
      INCLUDE 'params' 
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER NewEireneSurface_06

      INTEGER i1,i2,ik,ik1,ik2,side,sur1,sur2,sur3,iliin,ntmp,
     .        type,index1,index2,ilside,ilswch,region,code
      REAL    x1,x2,xcen,y1,y2,z1,z2,ycen,angle,dangle,rad,ewall,
     .        material,recycf,recyct

      nsurface = 0

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
          surface(nsurface)%index(1) = i1               ! Index of VESSEL_WALL surface in the fluid code wall array
          surface(nsurface)%v(1,1) =  wallpt(i2,20)     ! x coordinate of side vertex 1 (m)
          surface(nsurface)%v(2,1) =  wallpt(i2,21)     ! y coordinate
          surface(nsurface)%v(3,1) = -1.0E+18           ! z (code currently assumes toroidal symmetry)
          surface(nsurface)%v(1,2) =  wallpt(i1,20)     ! x coordinate of side vertex 2 (m)
          surface(nsurface)%v(2,2) =  wallpt(i1,21)     ! y
          surface(nsurface)%v(3,2) =  1.0E+18           ! z
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

          surface(nsurface)%v(1,1) =  eirasdat(i1,3)
          surface(nsurface)%v(2,1) =  eirasdat(i1,4)
          surface(nsurface)%v(3,1) = -1.0E+18
          surface(nsurface)%v(1,2) =  eirasdat(i1+1,3)
          surface(nsurface)%v(2,2) =  eirasdat(i1+1,4)
          surface(nsurface)%v(3,2) =  1.0E+18

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
            surface(nsurface)%v(1,1) = x1
            surface(nsurface)%v(2,1) = y1
            surface(nsurface)%v(3,1) = eirasdat(i1,4)
            surface(nsurface)%v(1,2) = x2
            surface(nsurface)%v(2,2) = y2
            surface(nsurface)%v(3,2) = eirasdat(i1,5)
          ENDDO

        ELSEIF (eirasdat(i1,1).EQ.12.0) THEN
c...      Rectangular tube: 
          x1 = eirasdat(i1,2)
          y1 = eirasdat(i1,3)
          x2 = eirasdat(i1,4)
          y2 = eirasdat(i1,5)

          nsurface = NewEireneSurface_06(VESSEL_WALL)
          surface(nsurface)%index(2) = NINT(eirasdat(i1,10))  ! Additional surface index
          surface(nsurface)%v(1,1) = x1
          surface(nsurface)%v(2,1) = y1
          surface(nsurface)%v(3,1) = eirasdat(i1,6)
          surface(nsurface)%v(1,2) = x2
          surface(nsurface)%v(2,2) = y1
          surface(nsurface)%v(3,2) = eirasdat(i1,7)

          nsurface = NewEireneSurface_06(VESSEL_WALL)
          surface(nsurface)%index(2) = NINT(eirasdat(i1,10))  
          surface(nsurface)%v(1,1) = x2
          surface(nsurface)%v(2,1) = y1
          surface(nsurface)%v(3,1) = eirasdat(i1,6)
          surface(nsurface)%v(1,2) = x2
          surface(nsurface)%v(2,2) = y2
          surface(nsurface)%v(3,2) = eirasdat(i1,7)

          nsurface = NewEireneSurface_06(VESSEL_WALL)
          surface(nsurface)%index(2) = NINT(eirasdat(i1,10))  
          surface(nsurface)%v(1,1) = x2
          surface(nsurface)%v(2,1) = y2
          surface(nsurface)%v(3,1) = eirasdat(i1,6)
          surface(nsurface)%v(1,2) = x1
          surface(nsurface)%v(2,2) = y2
          surface(nsurface)%v(3,2) = eirasdat(i1,7)

          nsurface = NewEireneSurface_06(VESSEL_WALL)
          surface(nsurface)%index(2) = NINT(eirasdat(i1,10))  
          surface(nsurface)%v(1,1) = x1
          surface(nsurface)%v(2,1) = y2
          surface(nsurface)%v(3,1) = eirasdat(i1,6)
          surface(nsurface)%v(1,2) = x1
          surface(nsurface)%v(2,2) = y1
          surface(nsurface)%v(3,2) = eirasdat(i1,7)

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
            surface(nsurface)%v(1,1) =  x1
            surface(nsurface)%v(2,1) =  y1
            surface(nsurface)%v(3,1) =  z1
            surface(nsurface)%v(1,2) =  x2
            surface(nsurface)%v(2,2) =  y2
            surface(nsurface)%v(3,2) =  z2
          ENDDO
        
        ELSE
        ENDIF
      ENDDO

c...  Non-default standard surfaces on the magnetic grid:
c     *TEMP* move lower eventually -- here now to be compatible with block 3a surface index of core
c     boundary in Detlev's reference input file...
      nsurface = NewEireneSurface_06(NON_DEFAULT_STANDARD)
      surface(nsurface)%subtype  = MAGNETIC_GRID_BOUNDARY
      IF (cgridopt.EQ.LINEAR_GRID) THEN
        surface(nsurface)%surtxt   = '* core, reflecting (OSM)'
        surface(nsurface)%index(1) = 1         
        surface(nsurface)%index(2) = MAXNKS    
        surface(nsurface)%index(3) = 2         
        surface(nsurface)%index(4) = 14      
        surface(nsurface)%index(5) = 0
        surface(nsurface)%iliin  = 1         
        surface(nsurface)%ilcol  = 3         
      ELSE
        surface(nsurface)%surtxt   = '* core, absorbing (OSM)'
        surface(nsurface)%index(1) = 1       ! Cell index start location of surface
        surface(nsurface)%index(2) = MAXNKS  ! Cell index end   
        surface(nsurface)%index(3) = 2       ! Ring in the magnetic grid 
        surface(nsurface)%index(4) = 14      ! Radial surface of cell (1-4 cell side here)
        surface(nsurface)%index(5) = 0
        surface(nsurface)%iliin  = 2         ! Reflection type (ILIIN=2 is 100% absorbing)
        surface(nsurface)%ilcol  = 3         ! Colour of surface in EIRENE plots
      ENDIF

c...  Poloidal surfaces for Eirene strata (needs to be generalized to 
c     allow more than 2 surface strata):
      DO i1 = IKLO, IKHI 
        nsurface = NewEireneSurface_06(NON_DEFAULT_STANDARD)
        surface(nsurface)%subtype  = STRATUM
        surface(nsurface)%surtxt   = '* target (OSM)'
        IF (cgridopt.EQ.LINEAR_GRID) THEN
          surface(nsurface)%index(1) = irsep               ! Ring index start location of surface
          surface(nsurface)%index(2) = nrs-1               ! Ring index end
          surface(nsurface)%index(3) = i1                  ! Target (IKLO=inner, IKHI=outer)
          surface(nsurface)%index(5) = 0
        ELSE
          surface(nsurface)%index(1) = irsep               ! Ring index start location of surface
          surface(nsurface)%index(2) = nrs                 ! Ring index end
          surface(nsurface)%index(3) = i1                  ! Target (IKLO=inner, IKHI=outer)
          surface(nsurface)%index(5) = 0
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
      ENDDO

c...  (Core boundary surface should ideally be set here (and not above).)

c...  Setup SOL radial boundary surfaces:
      ik1 = 0
      ik2 = 0
      DO ik = 1, nks(irwall)-1
c...
        IF (ik1.EQ.0) ik1 = ik
        IF ((ikins(ik,irwall).NE.ikins(ik+1,irwall)-1).OR.
     .      (irins(ik,irwall).NE.irins(ik+1,irwall)  )) ik2 = ik
        IF (ik.EQ.nks(irwall)-1.AND.ik2.EQ.0) ik2 = nks(irwall)

        IF (ik2.NE.0) THEN
          IF (irouts(ikins(ik1,irwall),irins(ik1,irwall)).EQ.
     .        irwall) THEN
            side = 23
          ELSE
            side = 14
          ENDIF
          nsurface = NewEireneSurface_06(NON_DEFAULT_STANDARD)
          surface(nsurface)%subtype  = MAGNETIC_GRID_BOUNDARY
          surface(nsurface)%surtxt   = '* radial SOL (OSM)'
          surface(nsurface)%index(1) = ikins(ik1,irwall)  ! Cell index start location of surface       
          surface(nsurface)%index(2) = ikins(ik2,irwall)  ! Cell index end                             
          surface(nsurface)%index(3) = irins(ik1,irwall)  ! Ring in the magnetic grid                  
          surface(nsurface)%index(4) = side               ! Radial surface of cell
          surface(nsurface)%index(5) = 0                                                                
          surface(nsurface)%iliin  = -1                   ! Reflection type (ILIIN=-1 is 100% transparent)
          surface(nsurface)%ilcol  = 4                    ! Colour of surface in EIRENE plots
          ik1 = 0
          ik2 = 0
        ENDIF
      ENDDO

c...  PFZ radial boundary:  
      IF (cgridopt.EQ.LINEAR_GRID.OR.irtrap.GT.nrs) THEN
      ELSE
        nsurface = NewEireneSurface_06(NON_DEFAULT_STANDARD)
        surface(nsurface)%subtype  = MAGNETIC_GRID_BOUNDARY
        surface(nsurface)%surtxt   = '* radial PFR (OSM)'
        surface(nsurface)%index(1) = 1                         ! Cell index start location of surface       
        surface(nsurface)%index(2) = MAXNKS                    ! Cell index end                             
        surface(nsurface)%index(3) = irtrap + 1                ! Ring in the magnetic grid                  
        surface(nsurface)%index(4) = 14                        ! Radial surface of cell (1-4 cell side here)
        surface(nsurface)%index(5) = 0                                                                      
        surface(nsurface)%iliin  = -1                          ! Reflection type (ILIIN=-1 is 100% transparent)
        surface(nsurface)%ilcol  = 4                           ! Colour of surface in EIRENE plots
      ENDIF

c...  Over-ride default surface properties:
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

c          WRITE(0,*) '-->',surface(i1)%index(1),surface(i1)%index(2),
c     .      surface(i1)%type,index1,index2

          IF ((surface(i1)%type   .EQ.NON_DEFAULT_STANDARD  .AND.   ! Code not tested as I forgot that
     .         surface(i1)%subtype.EQ.MAGNETIC_GRID_BOUNDARY.AND.   ! I was using the template directly
     .           type.EQ.1.AND.                                 ! for the case in question... *sigh*
     .           surface(i1)%index(3).GE.index1.AND.
     .           surface(i1)%index(3).LE.index2).OR.  ! *whew!*
     .        (surface(i1)%type.EQ.VESSEL_WALL.AND.
     .         ((type.EQ.2.AND.
     .           surface(i1)%index(1).GE.index1.AND.
     .           surface(i1)%index(1).LE.index2).OR.
     .          (type.EQ.3.AND.
     .           surface(i1)%index(2).GE.index1.AND.
     .           surface(i1)%index(2).LE.index2)))) THEN

            surface(i1)%iliin  = NINT(eirspdat(i2,3))
            surface(i1)%ilside = NINT(eirspdat(i2,4))
            surface(i1)%ilswch = NINT(eirspdat(i2,5))
            surface(i1)%recyct = eirspdat(i2,8)

c            WRITE(0,*) 'IL:',i1,surface(i1)%type
c            WRITE(0,*) '   ',surface(i1)%iliin,surface(i1)%ilside,
c     .                       surface(i1)%ilswch
          ENDIF
        ENDDO
      ENDDO
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
      DO i1 = 1, ntmp
        IF (surface(i1)%type.NE.VESSEL_WALL) CYCLE
c...    Surface matching criteria:
        iliin    = surface(i1)%iliin
        ilside   = surface(i1)%ilside
        ilswch   = surface(i1)%ilswch
        material = surface(i1)%material
        ewall    = surface(i1)%ewall
        recycf   = surface(i1)%recycf
        recyct   = surface(i1)%recyct
c...    Check if a non-default standard surface has already been defined
c       that matches the above citeria:
        sur3 = 0
        DO i2 = ntmp+1, nsurface
          IF (surface(i2)%iliin   .EQ.iliin   .AND.
     .        surface(i2)%ilside  .EQ.ilside  .AND.
     .        surface(i2)%ilswch  .EQ.ilswch  .AND.
     .        surface(i2)%material.EQ.material.AND.
     .        surface(i2)%ewall   .EQ.ewall   .AND.
     .        surface(i2)%recycf  .EQ.recycf  .AND.
     .        surface(i2)%recyct  .EQ.recyct) sur3 = i2
        ENDDO
        IF (sur3.EQ.0) THEN
c...      Existing surface matching the criteria for the current surface 
c         not found, so add a surface:
          sur2 = sur2 + 1
          nsurface = NewEireneSurface_06(NON_DEFAULT_STANDARD)
          surface(nsurface)%subtype  = ADDITIONAL
          surface(nsurface)%iliin    = surface(i1)%iliin
          surface(nsurface)%ilside   = surface(i1)%ilside
          surface(nsurface)%ilswch   = surface(i1)%ilswch
          surface(nsurface)%index(5) = sur2
          surface(nsurface)%reflect  = LOCAL
          surface(nsurface)%material = surface(i1)%material
          surface(nsurface)%ewall    = surface(i1)%ewall
          surface(nsurface)%recyct   = surface(i1)%recyct
          surface(nsurface)%recycf   = surface(i1)%recycf
          SELECTCASE (iliin)
            CASE (0)
              surface(nsurface)%surtxt = 
     .          '* transparent non-switching surface (DIVIMP)'
            CASE (1)
              surface(nsurface)%surtxt = '* material surface (DIVIMP)'      
            CASE (2)                
              surface(nsurface)%surtxt = '* pumping surface (DIVIMP)'      
            CASE DEFAULT
              CALL ER('DefineEireneSurfaces','Invalid ILIIN',*99)
          ENDSELECT
c...      Assign the wall surface to the non-default standard surface:
          surface(i1)%index(3) = sur2
        ELSE
c...      Assign the current wall surface to the identified non-default
c         standard surface that already exists:
          surface(i1)%index(3) = surface(sur3)%index(5) 
        ENDIF
      ENDDO

c...  Assign block 3a surface index to non-default standard surfaces:
      i2 = 0
      DO i1 = 1, nsurface
        IF (surface(i1)%type.EQ.NON_DEFAULT_STANDARD) THEN
          i2 = i2 + 1
          surface(i1)%num = i2
        ENDIF
      ENDDO

      RETURN
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
      USE mod_eirene06
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'


      REAL       TOL
      PARAMETER (TOL=1.0E-07)

      REAL GetMach,GetJsat,GetFlux 

      INTEGER ike,ik,ir,id,i1,in,it,region
      REAL    fact,Bfrac,e_pot(0:MAXNKS+1,MAXNRS),frac,maxepot
      REAL*8  Bx,By,Bz,beta,brat,deltax,deltay,x(3),y(3)


      WRITE(0,*) 'PROCESSING MAGNETIC GRID'


c...  Rough crack at e-potential:
      e_pot = 0.0
      DO ir = irsep, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        e_pot(0,ir) = 0.0
        DO ik = 1, nks(ir)
          IF (ik.EQ.1) THEN
            e_pot(ik,ir) = -0.5 * kes(ik,ir) * kss(ik,ir)
          ELSE
            e_pot(ik,ir) = e_pot(ik-1,ir) -
     .                     kes(ik,ir) * (kss(ik,ir) - kss(ik-1,ir))
          ENDIF
          IF (ik.EQ.nks(ir)) e_pot(ik+1,ir) = e_pot(ik,ir) -
     .                         kes(ik,ir) * (ksmaxs(ir) - kss(ik,ir))
        ENDDO
c        DO ik = 0, nks(ir)+1
c          WRITE(0,*) 'E_POT:',ik,ir,e_pot(ik,ir)
c        ENDDO  
        DO ik = 1, nks(ir)
          frac = kss(ik,ir) / ksmaxs(ir)
          e_pot(ik,ir) = e_pot(ik,ir) - frac * e_pot(nks(ir)+1,ir)
        ENDDO  
        e_pot(nks(ir)+1,ir) = 0.0
c        DO ik = 0, nks(ir)+1
c          WRITE(0,*) 'E_POT:',ik,ir,e_pot(ik,ir)
c        ENDDO  
c        STOP 'sdfsd'
      ENDDO
      ir = irsep
      maxepot = -1.0E+20
      DO ik = 1, nks(ir)
        IF (e_pot(ik,ir).GT.maxepot) maxepot = e_pot(ik,ir)
      ENDDO
      DO ir = 2, irsep-1
        e_pot(1:nks(ir),ir) = maxepot
      ENDDO

c...  Load cell geometry and volume quantities for the magnetic/fluid grid:
c      fact = qtim * qtim * emi / crmi
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
          cell(ncell)%type = 1
c...      Cell indices on the magnetic grid:
          cell(ncell)%index(1) = ik                               ! Cell index
          cell(ncell)%index(2) = ir                               ! Ring index
c...      Radial cell surfaces:
          cell(ncell)%sideindex(1,2) = 23                         ! Side index for cell surface 2
          cell(ncell)%sideindex(1,4) = 14                         ! Side index for cell surface 4
c...      Poloidal surfaces of note (targets):
          cell(ncell)%sideindex(2,1:4) = 0
          IF (ir.GE.irsep) THEN
            IF (ik.EQ.1      ) cell(ncell)%sideindex(2,1) = IKLO  ! Target index for cell surface 1
            IF (ik.EQ.nks(ir)) cell(ncell)%sideindex(2,3) = IKHI  ! Target index for cell surface 3
          ENDIF
c...      Cell vertices (4 vertices assumed at present):            
          id = korpg(ik,ir)
          DO i1 = 1, nvertp(id)
            cell(ncell)%r(i1) = DBLE(rvertp(i1,id))               ! x (m)
            cell(ncell)%z(i1) = DBLE(zvertp(i1,id))               ! y
          ENDDO
c...      B-field components (approximate):         
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
          By = brat * DSQRT(1.0D0/(1.0D0+beta**2)) * DSIGN(1.0D0,deltay)
          Bx = beta * By
c...      CBPHI is the on-axis B-field value specified in the OSM input file:
          Bfrac = cbphi * r0 / rs(ik,ir)                         ! Rough scaling of B-field
c...      Plasma quantities:
          cell(ncell)%plasma(1) = ktebs(ik,ir)                   ! Te (eV)
          cell(ncell)%plasma(2) = ktibs(ik,ir)                   ! Ti (eV)
          cell(ncell)%plasma(3) = knbs (ik,ir)                   ! ni (eV) (ne=ni assumed at present)
          cell(ncell)%plasma(4) = SNGL(Bx) * kvhs(ik,ir)         ! vx (m-1 s-1)
          cell(ncell)%plasma(5) = SNGL(By) * kvhs(ik,ir)         ! vy
          cell(ncell)%plasma(6) = SNGL(Bz) * kvhs(ik,ir)         ! vz
c...      E&M quantities:
          cell(ncell)%bfield(1) = SNGL(Bx) * Bfrac               ! Bx (Tesla) (normalized on Eirene side)
          cell(ncell)%bfield(2) = SNGL(By) * Bfrac               ! By 
          cell(ncell)%bfield(3) = SNGL(Bz) * Bfrac               ! Bz 
          cell(ncell)%efield(1) = SNGL(Bx) * kes(ik,ir) / fact   ! Ex (not required by EIRENE)
          cell(ncell)%efield(2) = SNGL(By) * kes(ik,ir) / fact   ! Ey
          cell(ncell)%efield(3) = SNGL(Bz) * kes(ik,ir) / fact   ! Ez

          cell(ncell)%e_pot = e_pot(ik,ir)                       ! Electric potential (estimate)
        ENDDO
      ENDDO


      IF (tetrahedrons) THEN
c   *** HACK ***
        cell(1)%plasma(1) = 10.0                  ! Te (eV)
        cell(1)%plasma(2) = 10.0                  ! Ti (eV)
        cell(1)%plasma(3) = 5.0E+21               ! ni (eV) (ne=ni assumed at present)
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
        ENDDO
      ENDDO      
      ntardat = it

      WRITE(0,*) 'DONE'

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
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      INTEGER i1,type,nspez,nasor,i2,insor,istrata
      LOGICAL assign_LO,assign_HI,assign_volrec

      TYPE(type_strata) :: tmpstrata      

      WRITE(0,*) 'NSTRATA:',nstrata,strata(1)%type
c      STOP
c      nstrata = 0

c...  Decide if default strata should be assigned:
      assign_LO     = .TRUE.
      assign_HI     = .TRUE.
      assign_volrec = .TRUE.
      DO istrata = 1, nstrata 
        IF (strata(istrata)%type         .EQ.1.0.AND.
     .      strata(istrata)%range_cell(1).EQ.-1 ) assign_LO     =.FALSE.  ! -1 = LO target, -2 = high target
        IF (strata(istrata)%type         .EQ.1.0.AND.
     .      strata(istrata)%range_cell(1).EQ.-2 ) assign_HI     =.FALSE.
        IF (strata(istrata)%type         .EQ.2.0) assign_volrec =.FALSE.
      ENDDO

      WRITE(0,*) 'STRATA:',assign_LO,assign_HI,assign_volrec

c...  Low IK target:
      IF (assign_LO) THEN
        nstrata = nstrata + 1

        strata(nstrata)%range_cell(1) = -1  ! No idea what this is for, also there for ASSIGN_HI...

        strata(nstrata)%type    = 1.0
        strata(nstrata)%indsrc  = 1
        strata(nstrata)%txtsou  = '* D+ bulk ions, low index target'
c        strata(nstrata)%npts    = 100
        strata(nstrata)%npts    = -90000
        strata(nstrata)%ninitl  = -1
        strata(nstrata)%nemods  =  3
        strata(nstrata)%flux    = 1.0
        strata(nstrata)%species_tag = 'FFFT'
        strata(nstrata)%nspez   =  1
        strata(nstrata)%distrib = 'FFTFF'
        strata(nstrata)%inum    =  1
        strata(nstrata)%indim   =  4
        strata(nstrata)%insor   = -2
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
        strata(nstrata)%ninitl  = -1
        strata(nstrata)%nemods  =  3
        strata(nstrata)%flux    = 1.0
        strata(nstrata)%species_tag = 'FFFT'
        strata(nstrata)%nspez   =  1
        strata(nstrata)%distrib = 'FFTFF'
        strata(nstrata)%inum    =  1
        strata(nstrata)%indim   =  4
        strata(nstrata)%insor   = -3
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
        strata(nstrata)%npts    = 100
c        strata(nstrata)%npts    = -90000
        strata(nstrata)%ninitl  = -1
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
      DO istrata = 1, nstrata
        WRITE(0,*) 'STRATA:TYPE=',NINT(strata(istrata)%type)

        SELECTCASE (NINT(strata(istrata)%type))
          CASE (999)  ! For compatibility with old strata definition format, see below
          CASE (1)
          CASE (2)
          CASE (3)
            IF (strata(istrata)%type.EQ.3.1) THEN  !***(small) hack***
              nspez = 2
              insor = 2 ! 1
              nasor = 1
              beam  = 1 ! Turn on beams in input file...
            ELSE
              nspez = 1
              insor = 1
              nasor = 0
            ENDIF
            strata(istrata)%indsrc  = -1
            strata(istrata)%txtsou  = '* point injection, '//
     .                                strata(istrata)%txtsou
            strata(istrata)%ninitl  = -1
            strata(istrata)%nemods  =  1
            strata(istrata)%flux    = strata(istrata)%flux *
     .                                strata(istrata)%flux_fraction
            strata(istrata)%species_tag = 'FFFF'
            i2 = strata(istrata)%species
            strata(istrata)%species_tag(i2:i2) = 'T'
            strata(istrata)%nspez   = nspez
            strata(istrata)%distrib = 'TFFFF'
            strata(istrata)%inum    = 0
            strata(istrata)%indim   = 0
            strata(istrata)%insor   = insor
            strata(istrata)%sorwgt  = 1.0
            strata(istrata)%sorlim  = 220.0
            strata(istrata)%sorind  = 0.0
            strata(istrata)%nrsor   = 0
            strata(istrata)%nasor   = nasor
            strata(istrata)%sorad(1:3) = strata(istrata)%sorad(1:3) *  ! Convert to cm
     .                                   100.0  

c            WRITE(0,*) strata(istrata)%sorad(1:3)
c            STOP 'sdgsdgd'
            strata(istrata)%soreni  = 0.0
          CASE (4)  ! Puff wall surface (not working...)
            WRITE(0,*) 'WALL SURFACE NEUTRAL INJECTION NOT WORKING'
            STOP 
          CASE DEFAULT
            CALL ER('SetupEireneStrata','Unrecognized strata type',*99)
        ENDSELECT

      ENDDO





c...  OLD SPECIFICATION: User specified neutral injection/puffing:
      DO i1 = 1, eirnpuff

c        STOP 'OLD PUFF SPECIFICATION DEACTIVATED'

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
            strata(nstrata)%ninitl  = -1
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
            strata(nstrata)%ninitl  = -1                      
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
      DO istrata = 1, nstrata
        IF (strata(istrata)%type.GE.3.0) THEN
          tmpstrata = strata(istrata)
          DO i1 = istrata, nstrata-1
            strata(i1) = strata(i1+1)
          ENDDO
          strata(nstrata) = tmpstrata  
        ENDIF
      ENDDO


      WRITE(0,*) 'STRATA:',nstrata,strata(1:nstrata)%type
      WRITE(0,*) 'STRATA:',strata(1:nstrata)%type


c      STOP 'sdgsgsdsdgsgd'

c...  Debugging:
c      strata(1:3)%npts = 2
c      strata(1:10)%ninitl = 11111 ! 22222 ! 99887

      RETURN
 99   WRITE(0,*) '  INDEX =',istrata
      WRITE(0,*) '  TYPE  =',strata(istrata)%type
      STOP
      END
c
c ======================================================================
c ======================================================================
c
      SUBROUTINE ReadEireneResults_06
      USE mod_eirene06
      IMPLICIT none

      CALL ReadParticleTracks_04

      CALL LoadEireneData_06
     
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
      SUBROUTINE LoadEireneData_06
      USE mod_geometry
      USE mod_eirene06_parameters
      USE mod_eirene06
      USE mod_eirene06_locals
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER fp,ntally,ndata,icount,index(30),ik,ir,i1,
     .        iblk,iatm,imol,iion,ipho,ilin,isur
      LOGICAL goodeof,output
      REAL    rdum(30),frac,norm,
     .        sumion,amps,pflux
      CHARACTER buffer*256,species*32

      INTEGER iobj,igrp
      INTEGER*2, ALLOCATABLE :: fluid_ik(:),fluid_ir(:)
      REAL     , ALLOCATABLE :: tdata(:,:,:)

      output = .FALSE.

      ALLOCATE(tdata(MAXNKS,MAXNRS,5))
      tdata = 0.0

      IF (tetrahedrons) THEN
        ALLOCATE(fluid_ik(nobj))
        ALLOCATE(fluid_ir(nobj))
        fluid_ik = 0
        fluid_ir = 0
        DO iobj = 1, nobj
          igrp = obj(iobj)%group
          IF (grp(igrp)%origin.EQ.GRP_MAGNETIC_GRID) THEN
            fluid_ik(iobj) = obj(iobj)%index(IND_IK)
            fluid_ir(iobj) = obj(iobj)%index(IND_IR)
          ENDIF
        ENDDO
      ELSE
        ALLOCATE(fluid_ik(ntri))
        ALLOCATE(fluid_ir(ntri))
        fluid_ik = 0
        fluid_ir = 0
        DO iobj = 1, ntri
          IF (tri(iobj)%type.EQ.MAGNETIC_GRID) THEN  
            fluid_ik(iobj) = tri(iobj)%index(1)                   ! Should pull these from .transfer
            fluid_ir(iobj) = tri(iobj)%index(2)   
          ENDIF
        ENDDO
      ENDIF

      goodeof = .FALSE.

c      pinion = 0.0   ! TEMP... 
c      pinrec = 0.0

      pinalpha = 0.0
      pinline = 0.0
      pinatom = 0.0
      pinmol = 0.0

      hescpd = 0.0
      hescal = 0.0

      IF (rel_opt.EQ.1.OR.rel_opt.EQ.3) THEN
        frac = rel_frac
      ELSE
        frac = 1.0
      ENDIF
      IF (sloutput)
     .  WRITE(0     ,*) 'RELAXATION FRACTION FOR EIRENE06:',frac
      WRITE(PINOUT,*) 'RELAXATION FRACTION FOR EIRENE06:',frac

      fp = 99
      OPEN(UNIT=fp,FILE='eirene.transfer',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)


      iblk = 0
      iatm = 0
      imol = 0
      ipho = 0
      ilin = 0
      DO WHILE (.TRUE.)
        READ(fp,'(A256)',END=10) buffer
        IF     (buffer(1:16).EQ.'* BULK PARTICLES') THEN
          IF (output) WRITE(0,*) '===BULK PARTICLES==='
          iblk = iblk + 1

          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         ! Check...
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          

          WRITE(0,*) 'NDATA:',ndata,nobj

c          READ(fp,'(A256)',END=10) buffer
c          WRITE(0,*) 'BULK1:',TRIM(buffer)
c          READ(buffer,*,ERR=97) ntally
c          READ(fp,'(A256)',END=10) buffer
c          WRITE(0,*) 'BULK2:',TRIM(buffer)
c          READ(buffer,*,ERR=97) ndata                         
c          READ(fp,'(A256)',END=10) buffer
c          WRITE(0,*) 'BULK3:',TRIM(buffer)
c          READ(buffer,*,ERR=97) (index(i1),i1=1,ntally)          

          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (fluid_ik(icount).NE.0) THEN
              ik = fluid_ik(icount)                          ! Should pull these from .transfer
              ir = fluid_ir(icount)
              IF (iblk.EQ.1) THEN                           ! Data for D+ only:
                tdata(ik,ir,1) = tdata(ik,ir,1) + rdum(5)   ! PINION
                tdata(ik,ir,2) = tdata(ik,ir,2) + rdum(6)   ! PINREC (relax?)
                tdata(ik,ir,3) = tdata(ik,ir,3) + rdum(7)   ! PINMP
                tdata(ik,ir,4) = tdata(ik,ir,4) + rdum(8)   ! PINQi
                tdata(ik,ir,5) = tdata(ik,ir,5) + rdum(9)   ! PINQe
              ENDIF
            ENDIF
          ENDDO
          IF (output) WRITE(0,*) '===DONE==='
        ELSEIF (buffer(1:12).EQ.'* TEST ATOMS') THEN
          IF (output) WRITE(0,*) '===TEST ATOMS==='
          iatm = iatm + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (fluid_ik(icount).NE.0) THEN
              ik = fluid_ik(icount)                          ! Should pull these from .transfer
              ir = fluid_ir(icount)
              IF (iatm.EQ.1) THEN                            ! Only for D, presumably the 1st atom species, need check...
                pinatom(ik,ir) = pinatom(ik,ir) + rdum(1) 
              ENDIF
            ENDIF
          ENDDO
          IF (output) WRITE(0,*) '===DONE==='
        ELSEIF (buffer(1:16).EQ.'* TEST MOLECULES') THEN
          IF (output) WRITE(0,*) '===TEST MOLECULES==='
          imol = imol + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)
            IF (fluid_ik(icount).NE.0) THEN
              ik = fluid_ik(icount)                          ! Should pull these from .transfer
              ir = fluid_ir(icount)
              IF (imol.EQ.1) THEN                            ! Need check...
                pinmol(ik,ir) = pinmol(ik,ir) + rdum(1)
              ENDIF
            ENDIF
          ENDDO
          IF (output) WRITE(0,*) '===DONE==='
        ELSEIF (buffer(1:11).EQ.'* TEST IONS') THEN
        ELSEIF (buffer(1:14).EQ.'* TEST PHOTONS') THEN
        ELSEIF (buffer(1:14).EQ.'* LINE EMISSIO') THEN
          IF (output) WRITE(0,*) '===LINE EMISSION==='
          ilin = ilin + 1   
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata   
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0

          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum)

            IF (fluid_ik(icount).NE.0) THEN
              ik = fluid_ik(icount)                          ! Should pull these from .transfer
              ir = fluid_ir(icount)
              IF (ilin.EQ.1) THEN
                pinline(ik,ir,1:5,H_BALPHA)=pinline(ik,ir,1:5,H_BALPHA)+ 
     .                                      rdum(1:5)
                pinline(ik,ir,6  ,H_BALPHA)=pinline(ik,ir,6  ,H_BALPHA)+ 
     .                                      rdum(7)
              ELSEIF (ilin.EQ.2) THEN
                pinline(ik,ir,1:6,H_BGAMMA)=pinline(ik,ir,1:6,H_BGAMMA)+ 
     .                                      rdum(1:6)
              ENDIF  
            ENDIF
          ENDDO
          IF (output) WRITE(0,*) '===DONE==='
        ELSEIF (buffer(1:6 ).EQ.'* MISC') THEN
c...      Check volumes:
        ELSEIF (buffer(1:13).EQ.'* PUMPED FLUX') THEN
          IF (output) WRITE(0,*) '===PUMPED FLUX==='
          DO WHILE (.TRUE.) 
            READ(fp,'(A256)',END=97,ERR=97) buffer
            IF (buffer(1:6).EQ.'* DONE') THEN
              goodeof = .TRUE.
              EXIT
            ENDIF
            READ(buffer,*) isur,species,amps
            IF     (species(1:6).EQ.'D     '.OR.
     .              species(1:6).EQ.'D(N=1)') THEN
              pflux = amps / ECH 
              hescpd  = hescpd + pflux
              IF (isur.NE.-1) hescal = hescal + pflux       ! Core ring assumption!!!!
            ELSEIF (species(1:6).EQ.'D2    ') THEN
              pflux = amps / ECH * 2.0
              hescpd  = hescpd + pflux
              IF (isur.NE.-1) hescal = hescal + pflux
            ELSE
              CALL WN('LoadEireneData','Unknown particle type')
            ENDIF          
          ENDDO
          IF (output) WRITE(0,*) '===DONE==='
        ELSEIF (buffer(1:1 ).EQ.'*') THEN
        ELSE
        ENDIF

      ENDDO
 10   CONTINUE

      CLOSE (fp)

c...  Relaxation into OSM arrays only?

c...  Normalize quantities (need to be careful vis-a-vis relaxation):
      sumion = 0.0
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE

        DO ik = 1, nks(ir)
          norm = 1.0 / kvols(ik,ir)
c          norm = 1.0 / (kvols(ik,ir) * eirtorfrac)

c...      Volume normalization:
          tdata(ik,ir,1:5) = tdata(ik,ir,1:5) / kvols(ik,ir)

c...      Linear relaxation:
          pinion(ik,ir) = (1.0-frac)*pinion(ik,ir) + frac*tdata(ik,ir,1)
          pinrec(ik,ir) = (1.0-frac)*pinrec(ik,ir) + frac*tdata(ik,ir,2)
          pinmp (ik,ir) = (1.0-frac)*pinmp (ik,ir) + frac*tdata(ik,ir,3)
          pinqi (ik,ir) = (1.0-frac)*pinqi (ik,ir) + frac*tdata(ik,ir,4)
          pinqe (ik,ir) = (1.0-frac)*pinqe (ik,ir) + frac*tdata(ik,ir,5)

c          pinion (ik,ir) = (1.0 - frac) * pinion(ik,ir) + 
c     .                            frac  * tdata (ik,ir,1) / kvols(ik,ir)
c          pinrec (ik,ir) = (1.0 - frac) * pinrec(ik,ir) + 
c     .                                    tdata (ik,ir,2) / kvols(ik,ir)

c...      Not relaxed, volume normalize:
          pinatom(ik,ir) = pinatom(ik,ir) * norm  !/ kvols(ik,ir)
          pinmol (ik,ir) = pinmol (ik,ir) * norm  !/ kvols(ik,ir)
          pinline(ik,ir,1:6,H_BALPHA) = pinline(ik,ir,1:6,H_BALPHA) * !/
     .                                  norm  !kvols(ik,ir)
          pinline(ik,ir,1:6,H_BGAMMA) = pinline(ik,ir,1:6,H_BGAMMA) * !/
     .                                  norm  !kvols(ik,ir)

c          DO i1 = 1, 6
c            pinline(ik,ir,i1,H_BALPHA) = pinline(ik,ir,i1,H_BALPHA) /
c    .                                    kvols(ik,ir)
c            pinline(ik,ir,i1,H_BGAMMA) = pinline(ik,ir,i1,H_BGAMMA) /
c    .                                    kvols(ik,ir)
c          ENDDO
c...
          pinalpha(ik,ir) = pinline(ik,ir,6,H_BALPHA)

          sumion = sumion + pinion(ik,ir) * kvols(ik,ir) * eirtorfrac
        ENDDO
      ENDDO

      WRITE(0,*) 'SUMION:',sumion

      CALL OutputEIRENE(67,'WORKING WITH EIRENE06')

c...  Make sure that the whole EIRENE data file was there:
      IF (.NOT.goodeof)  CALL ER('LoadEireneData','Problem with '//
     .                           'eirene.transfer file',*99)

c...  Need to save data again since Dalpha has been loaded into OBJ:
      CALL SaveGeometryData('tetrahedrons.raw')

c...  Clear arrays: 
      DEALLOCATE(fluid_ik)
      DEALLOCATE(fluid_ir)
      DEALLOCATE(tdata)
      IF (ALLOCATED(tri)) DEALLOCATE(tri) ! Move earlier and pass index mapping to/from Eirene
      IF (ALLOCATED(obj)) DEALLOCATE(obj) ! *same* but need to make more selective
      IF (ALLOCATED(srf)) DEALLOCATE(srf) ! at some point since even fluid objects will be stored here... 
      IF (ALLOCATED(vtx)) DEALLOCATE(vtx)

      RETURN
 97   WRITE(0,*) 'ERROR: PROBLEM READING DATA TRANSFER FILE'
      STOP
 98   WRITE(0,*) 'WARNING: eirene.transfer DATA FILE NOT FOUND'
      RETURN
 99   STOP
      END
c
c ======================================================================
c

