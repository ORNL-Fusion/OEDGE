c
c ======================================================================
c
c  subroutine: WriteEireneFiles_04
c
c
c
c
c
      subroutine WriteEireneFiles_04
      USE mod_eirene04
      IMPLICIT none

      include 'params'
      include 'cgeom'
      include 'comtor'
      INCLUDE 'slcom'

      INTEGER ik,ir,in1,in2,i1,id,ik1
      LOGICAL savedtriangles
    
      DATA savedtriangles /.FALSE./
 
      SAVE


      CALL OutputData(86,'Before calling EIRENE')


c...  Some variables used in the routine to write the '04 files, which do not have
c     access to the standard DIVIMP common blocks:
c     (EDGE2D: required)

      fluid_code = 'OSM'

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

c     (EDGE2D: set OPACITY=0, PHOTONS=0, TRIM=1 if wanting to use TRIM 
c      data for fast ion/atom reflection, BGK=0)
      opacity = eiropacity
      photons = eirphoton
      trim    = eirtrim
      bgk     = eirbgk


c...  A couple of checks:
c     (EDGE2D: not required, OSM specific)
      IF (eirniter.EQ.0.AND.(photons.GT.1.OR.bgk.EQ.3))
     .  CALL ER('WriteEireneFiles','Bad EIRNITER',*99)

      IF (eiropacity.NE.0.AND.eirphoton.NE.0)
     .  WRITE(0,*) 'WARNING: EIROPACITY.NE.0 AND EIRPHOTON.NE.0'


c...  Main calls for generating lots of triangles:  
      IF (.TRUE.) THEN

        IF (savedtriangles) THEN
c         (EDGE2D: not required, but handy for loading a previously generated triangluar mesh)
          CALL LoadTriangles

          CALL ALLOC_CELL(MAXNKS*nrs,nrs*2)
          CALL ProcessFluidGrid
          CALL AssignPlasmaQuantities
          CALL DEALLOC_CELL

          CALL WriteTriangleFiles

        ELSE
c         (EDGE2D: required)
          WRITE(0,*) 'BUILDING TRIANGLES'
          CALL ALLOC_VERTEX  (10000)   ! NEED TO ADD ACTIVE BOUNDS CHECKING FOR ALL THESE!
          CALL ALLOC_SURFACE (300)
          CALL ALLOC_TRIANGLE(15000)

          CALL DefineEireneSurfaces

c         * The calls below can be reordered to make things a little more efficient *

          CALL ALLOC_CELL(MAXNKS*nrs,nrs*2)
c...      Loads cell geometry and plasma quantities from fluid grid:
          CALL ProcessFluidGrid
c...      Fills the voids with triangles:
          CALL BuildFluidGridTriangles
c...      Builds connection map, amongst other things:
          CALL ProcessTriangles_04(0)
c...      Sets triangle volume quantities from data stored in ProcessFluidGrid:
          CALL AssignPlasmaQuantities
          CALL DEALLOC_CELL

c...      Fills the voids with triangles:
          CALL WritePolyFile(eirntri,MAXNRS,eirtri)
c...      Builds connection map, amongst other things:
          CALL ProcessTriangles_04(0)

c...      Writes the .points, .sides, .map and .plasma files to be loaded
c         by EIRENE:
          CALL WriteTriangleFiles

c...      Dumps trianlges to a binary file, for use with LoadTriangles:
          CALL SaveTriangles
          savedtriangles = .TRUE.
          WRITE(0,*) 'DONE'
        ENDIF
 
        WRITE(0,*) 'WRITING THE EIRENE INPUT FILE'
        IF (eirdata.EQ.1) CALL WriteEireneInputFile_04
        WRITE(0,*) 'DONE'

      ELSE
c       (EDGE2D: not required, OSM specific)
        STOP 'OLD EIRENE COUPLING NO LONGER SUPPORTED'
        CALL DumpTriangles
        CALL WriteEireneInputFile_04
      ENDIF


c...  TMP: Calculate radial 'sources':
c     (EDGE2D: not required, OSM specific)
      IF     (s28cfpdrft.GE.1) THEN
        WRITE(0,*) 'CALCULATING RADIAL FLUX AFTER RETURNING FROM EIRENE'
        CALL CalcRadialDrift(-1)
        CALL CalcRadialDrift(-2)
      ELSEIF (s28cfpdrft.EQ.-1) THEN
c...    Turning it off:
        WRITE(0,*) 'TURING OFF RADIAL FLUX AFTER RETURNING FROM EIRENE'
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            osmcfpflx(ik,ir,2) = (1.0 - rel_frac) * osmcfpflx(ik,ir,2) 
          ENDDO
        ENDDO
        CALL CalcRadialDrift(-2)
      ENDIF

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: DefineEireneSurfaces
c
      SUBROUTINE DefineEireneSurfaces
      USE mod_eirene04
      IMPLICIT none

      INTEGER region,code

      INCLUDE 'params' 
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER NewEireneSurface

      INTEGER i1,i2,ik,ik1,ik2,side,sur1,sur2,sur3,iliin,ntmp,
     .        type,index1,index2
      REAL    x1,x2,xcen,y1,y2,z1,z2,ycen,angle,dangle,rad,ewall,
     .        material,recycf,recyct

c...  Initialize surface array list:
      nadd = 0

c...  Load vessel wall segments from the default DIVIMP vessel 
c     wall array (WALLPTS), clockwise wall is assumed:
c
c * Need to add that 
c
c     (EDGE2D: required)
      IF (.TRUE.) THEN
        DO i1 = wallpts, 1, -1
          IF (i1.EQ.wallpts) THEN
            i2 = 1
          ELSE
            i2 = i1 + 1
          ENDIF
          IF (wallpt(i1,18).EQ.0.0) THEN
            nadd = NewEireneSurface(VESSEL_WALL)
            add(nadd)%index(1) = i1               ! Index of VESSEL_WALL surface in the fluid code wall array
            add(nadd)%v(1,1) =  wallpt(i2,20)     ! x coordinate of side vertex 1 (m)
            add(nadd)%v(2,1) =  wallpt(i2,21)     ! y coordinate
            add(nadd)%v(3,1) = -1.0E+18           ! z (code currently assumes toroidal symmetry)
            add(nadd)%v(1,2) =  wallpt(i1,20)     ! x coordinate of side vertex 2 (m)
            add(nadd)%v(2,2) =  wallpt(i1,21)     ! y
            add(nadd)%v(3,2) =  1.0E+18           ! z
          ENDIF
        ENDDO
      ENDIF

c...  User specified additional surfaces (EIRASDAT array):
c     (EDGE2D: probably don't need something like this just yet...)
      IF (.TRUE.) THEN

        DO i1 = 1, eirnasdat
          IF (eirasdat(i1,1) .EQ.-1.0.OR.
c     .        eirasdat(i1,2) .NE. 1.0.OR.
c...          Only load surfaces with a specified additional surface index:
     .        eirasdat(i1,10).EQ. 0.0) CYCLE

c... * The 7.0 designation needs to change since not all the 7.0 info is being used *
          IF     (eirasdat(i1,1).EQ.7.0) THEN

            nadd = NewEireneSurface(VESSEL_WALL)

            add(nadd)%index(2) = NINT(eirasdat(i1,10))  ! Additional surface index

            add(nadd)%v(1,1) =  eirasdat(i1,3)
            add(nadd)%v(2,1) =  eirasdat(i1,4)
            add(nadd)%v(3,1) = -1.0E+18
            add(nadd)%v(1,2) =  eirasdat(i1+1,3)
            add(nadd)%v(2,2) =  eirasdat(i1+1,4)
            add(nadd)%v(3,2) =  1.0E+18

          ELSEIF (eirasdat(i1,1).EQ.11.0) THEN
c...        Tube (fun!):
            dangle = (eirasdat(i1,8) - eirasdat(i1,7)) / eirasdat(i1,9) 
            xcen = eirasdat(i1,2)
            ycen = eirasdat(i1,3)
            rad = eirasdat(i1,6)
            DO angle = eirasdat(i1,7), 0.9999*eirasdat(i1,8), dangle
              x1 = xcen + rad * COS(angle*PI/180.0)
              y1 = ycen + rad * SIN(angle*PI/180.0)
              x2 = xcen + rad * COS((angle+dangle)*PI/180.0)
              y2 = ycen + rad * SIN((angle+dangle)*PI/180.0)

              nadd = NewEireneSurface(VESSEL_WALL)
              add(nadd)%index(2) = NINT(eirasdat(i1,10))  ! Additional surface index

              add(nadd)%v(1,1) = x1
              add(nadd)%v(2,1) = y1
              add(nadd)%v(3,1) = eirasdat(i1,4)
              add(nadd)%v(1,2) = x2
              add(nadd)%v(2,2) = y2
              add(nadd)%v(3,2) = eirasdat(i1,5)
            ENDDO

          ELSEIF (eirasdat(i1,1).EQ.12.0) THEN
c...        Rectangular tube: 
            x1 = eirasdat(i1,2)
            y1 = eirasdat(i1,3)
            x2 = eirasdat(i1,4)
            y2 = eirasdat(i1,5)

            nadd = NewEireneSurface(VESSEL_WALL)
            add(nadd)%index(2) = NINT(eirasdat(i1,10))  ! Additional surface index

            add(nadd)%v(1,1) = x1
            add(nadd)%v(2,1) = y1
            add(nadd)%v(3,1) = eirasdat(i1,6)
            add(nadd)%v(1,2) = x2
            add(nadd)%v(2,2) = y1
            add(nadd)%v(3,2) = eirasdat(i1,7)

            nadd = NewEireneSurface(VESSEL_WALL)
            add(nadd)%index(2) = NINT(eirasdat(i1,10))  

            add(nadd)%v(1,1) = x2
            add(nadd)%v(2,1) = y1
            add(nadd)%v(3,1) = eirasdat(i1,6)
            add(nadd)%v(1,2) = x2
            add(nadd)%v(2,2) = y2
            add(nadd)%v(3,2) = eirasdat(i1,7)

            nadd = NewEireneSurface(VESSEL_WALL)
            add(nadd)%index(2) = NINT(eirasdat(i1,10))  

            add(nadd)%v(1,1) = x2
            add(nadd)%v(2,1) = y2
            add(nadd)%v(3,1) = eirasdat(i1,6)
            add(nadd)%v(1,2) = x1
            add(nadd)%v(2,2) = y2
            add(nadd)%v(3,2) = eirasdat(i1,7)

            nadd = NewEireneSurface(VESSEL_WALL)
            add(nadd)%index(2) = NINT(eirasdat(i1,10))  

            add(nadd)%v(1,1) = x1
            add(nadd)%v(2,1) = y2
            add(nadd)%v(3,1) = eirasdat(i1,6)
            add(nadd)%v(1,2) = x1
            add(nadd)%v(2,2) = y1
            add(nadd)%v(3,2) = eirasdat(i1,7)

          ELSEIF (eirasdat(i1,1).EQ.14.0) THEN

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

              nadd = NewEireneSurface(VESSEL_WALL)

              add(nadd)%index(2) = NINT(eirasdat(i1,10))  ! Additional surface index

              add(nadd)%v(1,1) =  x1
              add(nadd)%v(2,1) =  y1
              add(nadd)%v(3,1) =  z1
              add(nadd)%v(1,2) =  x2
              add(nadd)%v(2,2) =  y2
              add(nadd)%v(3,2) =  z2

            ENDDO
          
          ELSE
          ENDIF

        ENDDO

      ENDIF

c...  Non-default standard surfaces on the magnetic grid:
c     (EDGE2D: required)
      IF (.TRUE.) THEN

c       *TEMP* move lower eventually -- here now to be compatible with block 3a surface index of core
c       boundary in Detlev's reference input file...
        nadd = NewEireneSurface(NON_DEFAULT_STANDARD)
        add(nadd)%subtype  = MAGNETIC_GRID_BOUNDARY
        IF (cgridopt.EQ.LINEAR_GRID) THEN
          add(nadd)%surtxt   = '* core boundary, reflecting (OSM)'
          add(nadd)%index(1) = 1         
          add(nadd)%index(2) = MAXNKS    
          add(nadd)%index(3) = 2         
          add(nadd)%index(4) = 14      
          add(nadd)%index(5) = 0
          add(nadd)%iliin  = 1         
          add(nadd)%ilcol  = 3         
        ELSE
          add(nadd)%surtxt   = '* core boundary, absorbing (OSM)'
          add(nadd)%index(1) = 1       ! Cell index start location of surface
          add(nadd)%index(2) = MAXNKS  ! Cell index end   
          add(nadd)%index(3) = 2       ! Ring in the magnetic grid 
          add(nadd)%index(4) = 14      ! Radial surface of cell (1-4 cell side here)
          add(nadd)%index(5) = 0
          add(nadd)%iliin  = 2         ! Reflection type (ILIIN=2 is 100% absorbing)
          add(nadd)%ilcol  = 3         ! Colour of surface in EIRENE plots
        ENDIF

c...    Surfaces for strata (needs to be generalized to allow more 
c       than 2 surface strata):
        DO i1 = IKLO, IKHI 
          IF (.TRUE.) THEN
            nadd = NewEireneSurface(NON_DEFAULT_STANDARD)
            add(nadd)%subtype  = STRATUM
            add(nadd)%surtxt   = '* target (OSM)'
            IF (cgridopt.EQ.LINEAR_GRID) THEN
              add(nadd)%index(1) = irsep               ! Ring index start location of surface
              add(nadd)%index(2) = nrs-1               ! Ring index end
              add(nadd)%index(3) = i1                  ! Target (IKLO=inner, IKHI=outer)
              add(nadd)%index(5) = 0
            ELSE
              add(nadd)%index(1) = irsep               ! Ring index start location of surface
              add(nadd)%index(2) = nrs                 ! Ring index end
              add(nadd)%index(3) = i1                  ! Target (IKLO=inner, IKHI=outer)
              add(nadd)%index(5) = 0
            ENDIF

            add(nadd)%reflect = LOCAL                  ! Set surface reflection model to LOCAL
            add(nadd)%iliin  = 1
            add(nadd)%ilside = 0
            add(nadd)%ilswch = 0
            add(nadd)%iltor  = 0
            add(nadd)%ilcell = 0
            add(nadd)%ilcol  = 4
            add(nadd)%material = tmater                ! Set surface material
            add(nadd)%ewall = -ttemp * 1.38E-23 / ECH  ! Set temperature
          ENDIF
        ENDDO

      ENDIF


c...  Set radial boundary surfaces of the standard magnetic grid:
c     (EDGE2D: required)
      IF (.TRUE.) THEN
c...    Core boundary surface should ideally be set here (and not above):

c...    Radial SOL boundaries:
        ik1 = 0
        ik2 = 0
        DO ik = 1, nks(irwall)-1
          IF (ik1.EQ.0) ik1 = ik
          IF (irins(ik,irwall).NE.irins(ik+1,irwall)) ik2 = ik
          IF (ik.EQ.nks(irwall)-1.AND.ik2.EQ.0) ik2 = nks(irwall)
          IF (ik2.NE.0) THEN
            IF (irouts(ikins(ik1,irwall),irins(ik1,irwall)).EQ.
     .          irwall) THEN
              side = 23
            ELSE
              side = 14
            ENDIF

            nadd = NewEireneSurface(NON_DEFAULT_STANDARD)
            add(nadd)%subtype  = MAGNETIC_GRID_BOUNDARY
            add(nadd)%surtxt   = '* radial SOL (OSM)'
            add(nadd)%index(1) = ikins(ik1,irwall)  ! Cell index start location of surface       
            add(nadd)%index(2) = ikins(ik2,irwall)  ! Cell index end                             
            add(nadd)%index(3) = irins(ik1,irwall)  ! Ring in the magnetic grid                  
            add(nadd)%index(4) = side               ! Radial surface of cell
            add(nadd)%index(5) = 0                                                                
            add(nadd)%iliin  = -1                   ! Reflection type (ILIIN=-1 is 100% transparent)
            add(nadd)%ilcol  = 4                    ! Colour of surface in EIRENE plots

            ik1 = 0
            ik2 = 0
          ENDIF
        ENDDO

c...    PFZ boundary:  
        IF (cgridopt.EQ.LINEAR_GRID.OR.irtrap.GT.nrs) THEN
        ELSE
          nadd = NewEireneSurface(NON_DEFAULT_STANDARD)
          add(nadd)%subtype  = MAGNETIC_GRID_BOUNDARY
          add(nadd)%surtxt   = '* radial PFR (OSM)'
          add(nadd)%index(1) = 1                         ! Cell index start location of surface       
          add(nadd)%index(2) = MAXNKS                    ! Cell index end                             
          add(nadd)%index(3) = irtrap + 1                ! Ring in the magnetic grid                  
          add(nadd)%index(4) = 14                        ! Radial surface of cell (1-4 cell side here)
          add(nadd)%index(5) = 0                                                                      
          add(nadd)%iliin  = -1                          ! Reflection type (ILIIN=1 is 100% transparent)
          add(nadd)%ilcol  = 4                           ! Colour of surface in EIRENE plots
        ENDIF

      ENDIF





c...  *HACK* set pumping surface..
      DO i1 = 1, nadd
        IF (add(i1)%type    .EQ.VESSEL_WALL.AND.
     .      add(i1)%index(2).EQ.200.0      ) THEN
          add(i1)%iliin = 2
        ENDIF
      ENDDO


c...  Adjust surface properties:
      DO i1 = 1, nadd
        DO i2 = 1, eirnspdat
          type   = NINT(eirspdat(i2,1))
          index1 = NINT(eirspdat(i2,2))
          IF (eirspdat(MIN(i2+1,eirnspdat),1).EQ.0.0) THEN
            index2 = NINT(eirspdat(i2+1,2))
          ELSE
            index2 = index1
          ENDIF
          IF (.TRUE.) THEN

            IF (add(i1)%type.EQ.VESSEL_WALL.AND.
     .           ((type.EQ.2.AND.
     .             add(i1)%index(1).GE.index1.AND.
     .             add(i1)%index(1).LE.index2).OR.
     .            (type.EQ.3.AND.
     .             add(i1)%index(2).GE.index1.AND.
     .             add(i1)%index(2).LE.index2))) THEN

              add(i1)%iliin  = NINT(eirspdat(i2,3))
              add(i1)%ilside = NINT(eirspdat(i2,4))
              add(i1)%ilswch = NINT(eirspdat(i2,5))
              add(i1)%recyct = eirspdat(i2,8)

            ENDIF

          ENDIF
        ENDDO
      ENDDO

c...  Build non-default standard surfaces to represent vessel wall surfaces when
c     writing block 3 of the the EIRENE04 input file.  A different surface is
c     defined each time a different surface reflection model and/or surface 
c     property (material, temperature, etc.) is identified.  At present, only the
c     TYPE.EQ.NON_DEFAULT_STANDARD surfaces in the list of surfaces are used
c     when writing block 3a in the EIRENE input file (and block 3b currently has
c     no surfaces):
c
c     (EDGE2D: if the above surfaces are defined correctly, this piece of code
c              should be the same for all fluid codes)
      IF (.TRUE.) THEN

c...    Count the number of non-default standard surfaces already 
c       defined (magnetic grid boundaries specified above):
        sur1 = 0
        DO i1 = 1, nadd
          IF (add(i1)%type.EQ.NON_DEFAULT_STANDARD) sur1 = sur1 + 1
        ENDDO

        sur1 = sur1 + 1
        sur2 = sur1 - 1

        ntmp = nadd

        DO i1 = 1, ntmp
          IF (add(i1)%type.NE.VESSEL_WALL) CYCLE

c...      Surface matching criteria:
          iliin    = add(i1)%iliin
          material = add(i1)%material
          ewall    = add(i1)%ewall
          recycf   = add(i1)%recycf
          recyct   = add(i1)%recyct

c...      Check if a non-default standard surface has already been defined
c         that matches the above citeria:
          sur3 = 0
          DO i2 = ntmp+1, nadd
            IF (add(i2)%iliin   .EQ.iliin   .AND.
     .          add(i2)%material.EQ.material.AND.
     .          add(i2)%ewall   .EQ.ewall   .AND.
     .          add(i2)%recycf  .EQ.recycf  .AND.
     .          add(i2)%recyct  .EQ.recyct) sur3 = i2
          ENDDO

          IF (sur3.EQ.0) THEN
c...        Existing surface matching the criteria for the current surface 
c           not found, so add a surface:
            sur2 = sur2 + 1

            nadd = NewEireneSurface(NON_DEFAULT_STANDARD)
            add(nadd)%subtype  = ADDITIONAL
            add(nadd)%iliin    = add(i1)%iliin
            add(nadd)%index(5) = sur2

            add(nadd)%reflect = LOCAL
            add(nadd)%material = add(i1)%material
            add(nadd)%ewall    = add(i1)%ewall
            add(nadd)%recyct   = add(i1)%recyct
            add(nadd)%recycf   = add(i1)%recycf

            IF     (iliin.EQ.1) THEN
              add(nadd)%surtxt   = '* material surface (OSM)'      
            ELSEIF (iliin.EQ.2) THEN
              add(nadd)%surtxt   = '* pumping surface (OSM)'      
            ELSEIF (iliin.EQ.0) THEN
              add(nadd)%surtxt   = '* transparent non-switching '//
     .                             'surface (OSM)'
            ELSE
              CALL ER('DefineEireneSurfaces','Invalid ILIIN value',*99)
            ENDIF

c...        Assign the wall surface to the non-default standard surface:
            add(i1)%index(3) = sur2

          ELSE
c...        Assign the current wall surface to the identified non-default
c           standard surface that already exists:
            add(i1)%index(3) = add(sur3)%index(5) 
          ENDIF

        ENDDO

c...    Assign block 3a surface index to non-default standard surfaces:
        i2 = 0
        DO i1 = 1, nadd
          IF (add(i1)%type.EQ.NON_DEFAULT_STANDARD) THEN
            i2 = i2 + 1
            add(i1)%num = i2
          ENDIF
        ENDDO

      ENDIF


      RETURN
 99   STOP
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
      SUBROUTINE ProcessFluidGrid
      USE mod_eirene04
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'


      REAL       TOL
      PARAMETER (TOL=1.0E-07)

      REAL GetMach,GetJsat,GetFlux 

      INTEGER ike,ik,ir,id,i1,in,it,region
      REAL    fact,Bfrac
      REAL*8  Bx,By,Bz,beta,brat,deltax,deltay,x(3),y(3)


      WRITE(0,*) '  PROCESSING MAGNETIC GRID'


c      fact = qtim * qtim * emi / crmi
      fact = 1.0

      ncell = 0

c...  Load cell geometry and volume quantities for the magnetic grid:
c     (EDGE2D: required)
c     (NRS is the number of rings on the grid)
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE

c       (NKS(IR) is the number of cells on ring IR)
        IF (ir.LT.irsep) THEN
          ike = nks(ir) - 1
        ELSE
          ike = nks(ir)
        ENDIF

        DO ik = 1, ike
          ncell = ncell + 1
 
          cell(ncell)%type = 1
c...      Cell indeces on the magnetic grid:
          cell(ncell)%index(1) = ik                               ! Cell index
          cell(ncell)%index(2) = ir                               ! Rind index
c...      Radial cell surfaces:
          cell(ncell)%sideindex(1,2) = 23                         ! Side index for cell surface 2
          cell(ncell)%sideindex(1,4) = 14                         ! Side index for cell surface 4
c...      Poloidal surfaces of note (targets):
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
c...      B-field components:         
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
          Bz = DSQRT(1.0 - brat**2)                               ! * I think this algebra is ok, but check needed *
          By = brat * DSQRT(1.0D0/(1.0D0+beta**2)) * DSIGN(1.0D0,deltay)
          Bx = beta * By
c...      CBPHI is the on-axis B-field value specified in the OSM input file:
          Bfrac = cbphi * r0 / rs(ik,ir)                          ! Rough scaling of B-field
c...      Plasma quantities:
          cell(ncell)%plasma(1) = ktebs(ik,ir)                   ! Te (eV)
          cell(ncell)%plasma(2) = ktibs(ik,ir)                   ! Ti (eV)
          cell(ncell)%plasma(3) = knbs (ik,ir)                   ! ni (eV) (ne=ni assumed at present)
          cell(ncell)%plasma(4) = SNGL(Bx) * kvhs(ik,ir)         ! vx (m-1 s-1)
          cell(ncell)%plasma(5) = SNGL(By) * kvhs(ik,ir)         ! vy
          cell(ncell)%plasma(6) = SNGL(Bz) * kvhs(ik,ir)         ! vz
c...      E&M quantities:
          cell(ncell)%bfield(1) = SNGL(Bx) * Bfrac               ! Bx (Tesla) (normalized on EIRENE side)
          cell(ncell)%bfield(2) = SNGL(By) * Bfrac               ! By 
          cell(ncell)%bfield(3) = SNGL(Bz) * Bfrac               ! Bz 
          cell(ncell)%efield(1) = SNGL(Bx) * kes(ik,ir) / fact   ! Ex (not required by EIRENE)
          cell(ncell)%efield(2) = SNGL(By) * kes(ik,ir) / fact   ! Ey
          cell(ncell)%efield(3) = SNGL(Bz) * kes(ik,ir) / fact   ! Ez
        ENDDO
      ENDDO


c...  Need to load target quantities as well, but not done here yet:
c     (EDGE2D: required)
      it = 0
      DO region = IKLO, IKHI  !  (EDGE2D: IKLO=1, IKHI=2)

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
c          IF (cgridopt.EQ.LINEAR_GRID) THEN
c            WRITE(0,*) ' **** BOGUS TARGET FLUX TO EIRENE04 ****'
c            tardat(it,7)  = 1.0E+22 * ECH
c          ELSE
          tardat(it,7)  = GetFlux(region,ir) / eirtorfrac /eirsrcmul*ECH  ! ion flux to target for species 1 (Amps) (req.)
c          ENDIF
          tardat(it,8)  = ktids(in)                                       ! Ti (eV)  (required)
          tardat(it,9)  = knds(in)                                        ! ni (m-3) (required)
          tardat(it,10) = kvds(in)                                        ! v_para (m s-1)
          tardat(it,11) = GetMach(kvds(in),kteds(in),ktids(in))           ! Mach no.
          tardat(it,12) = GetJsat(kteds(in),ktids(in),knds(in),kvds(in))  ! jsat 

c...      Suppress target flux:
c          WRITE(0,*) 'SUPFLX:',region,ir,supflx(region,ir)
          IF (supflx(region,ir).EQ.1) tardat(it,7) =tardat(it,7)*1.0E-15
        ENDDO
      ENDDO      
      ntardat = it

      WRITE(0,*) '  NTARDAT=',ntardat
      WRITE(0,*) '  DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c ======================================================================
c
      SUBROUTINE ReadEireneResults_04
      USE mod_eirene04
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'pindata'

      CALL ReadParticleTracks_04

      CALL LoadEireneData_04
     












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


c...  Assign FLXHWx arrays:
c
c     FLUXHW - FLUX OF HYDROGEN (ATOMS AND MOLECULES) TO THE WALL
c     FLXHW2 - FLUX OF HYDROGEN (ATOMS AND IONS) TO THE WALL
c     FLXHW3 - FLUX OF IMPURITIES SPUTTERED FROM THE WALL (N/A)
c     FLXHW4 - FLUX OF IMPURITIES REDEPOSITED ONTO THE WALL (N/A)
c     FLXHW5 - AVERAGE ENERGY OF ATOMS HITTING THE WALL (EV)
c     FLXHW6 - FLUX OF HYDROGEN ATOMS TO THE WALL
c     FLXHW7 - AVERAGE ENERGY OF MOLECULES HITTING THE WALL (eV)
c     FLXHW8 - EIRENE REPORTED HYDROGEN ION FLUXES TO THE WALL 
c
      CALL RZero(fluxhw,MAXSEG)
      CALL RZero(flxhw2,MAXSEG)
      CALL RZero(flxhw3,MAXSEG)
      CALL RZero(flxhw4,MAXSEG)
      CALL RZero(flxhw5,MAXSEG)
      CALL RZero(flxhw6,MAXSEG)
      CALL RZero(flxhw7,MAXSEG)
      CALL RZero(flxhw8,MAXSEG)



      RETURN
 99   STOP
      END
c
c ======================================================================
c ======================================================================
c
c
c ======================================================================
c
      SUBROUTINE LoadEireneData_04
      USE mod_eirene04
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER fp,ntally,ndata,icount,index(30),ik,ir,i1,
     .        iblk,iatm,imol,iion,ipho,ilin,isur
      LOGICAL goodeof,binary
      REAL    rdum(30),frac,
     .        sumion,amps,flux
      CHARACTER buffer*256,species*32

      REAL, ALLOCATABLE :: tdata(:,:,:)

      ALLOCATE(tdata(MAXNKS,MAXNRS,5))
      tdata = 0.0

      binary  = .FALSE.
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
     .  WRITE(0     ,*) 'RELAXATION FRACTION FOR EIRENE04:',frac
      WRITE(PINOUT,*) 'RELAXATION FRACTION FOR EIRENE04:',frac

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

          iblk = iblk + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         ! Check...
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum,binary)
            IF (tri(icount)%type.EQ.MAGNETIC_GRID) THEN  
              ik = tri(icount)%index(1)                   ! Should pull these from .transfer
              ir = tri(icount)%index(2)   
              IF (iblk.EQ.1) THEN                           ! Data for D+ only:
                tdata(ik,ir,1) = tdata(ik,ir,1) + rdum(5)   ! PINION
                tdata(ik,ir,2) = tdata(ik,ir,2) + rdum(6)   ! PINREC (relax?)
                tdata(ik,ir,3) = tdata(ik,ir,3) + rdum(7)   ! PINMP
                tdata(ik,ir,4) = tdata(ik,ir,4) + rdum(8)   ! PINQi
                tdata(ik,ir,5) = tdata(ik,ir,5) + rdum(9)   ! PINQe
              ENDIF
            ENDIF
          ENDDO

        ELSEIF (buffer(1:12).EQ.'* TEST ATOMS') THEN

          iatm = iatm + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum,binary)
            IF (tri(icount)%type.EQ.MAGNETIC_GRID) THEN  
              ik = tri(icount)%index(1)                   
              ir = tri(icount)%index(2)   
              IF (iatm.EQ.1) THEN                            ! Only for D, presumably the 1st atom species, need check...
                pinatom(ik,ir) = pinatom(ik,ir) + rdum(1) 
              ENDIF
            ENDIF
          ENDDO

        ELSEIF (buffer(1:16).EQ.'* TEST MOLECULES') THEN

          imol = imol + 1
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata                         
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum,binary)
            IF (tri(icount)%type.EQ.MAGNETIC_GRID) THEN  
              ik = tri(icount)%index(1)                   
              ir = tri(icount)%index(2)   
              IF (imol.EQ.1) THEN                            ! Need check...
                pinmol(ik,ir) = pinmol(ik,ir) + rdum(1)
              ENDIF
            ENDIF
          ENDDO

        ELSEIF (buffer(1:11).EQ.'* TEST IONS') THEN
        ELSEIF (buffer(1:14).EQ.'* TEST PHOTONS') THEN
        ELSEIF (buffer(1:14).EQ.'* LINE EMISSIO') THEN

          ilin = ilin + 1   
          READ(fp,*,ERR=97) ntally
          READ(fp,*,ERR=97) ndata   
          READ(fp,*,ERR=97) (index(i1),i1=1,ntally)          
          icount = 0
          DO WHILE (icount.LT.ndata)
            CALL NextLine(fp,ntally,icount,rdum,binary)
            IF (tri(icount)%type.EQ.MAGNETIC_GRID) THEN  
              ik = tri(icount)%index(1)                  
              ir = tri(icount)%index(2)   
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

        ELSEIF (buffer(1:6 ).EQ.'* MISC') THEN
c...      Check volumes:
        ELSEIF (buffer(1:13).EQ.'* PUMPED FLUX') THEN

          DO WHILE (.TRUE.) 
            READ(fp,'(A256)',END=97,ERR=97) buffer
            IF (buffer(1:6).EQ.'* DONE') THEN
              goodeof = .TRUE.
              EXIT
            ENDIF
            READ(buffer,*) isur,species,amps
            IF     (species(1:6).EQ.'D     '.OR.
     .              species(1:6).EQ.'D(N=1)') THEN
              flux = amps / ECH 
              hescpd  = hescpd + flux
              IF (isur.NE.-1) hescal = hescal + flux       ! Core ring assumption!!!!
            ELSEIF (species(1:6).EQ.'D2    ') THEN
              flux = amps / ECH * 2.0
              hescpd  = hescpd + flux
              IF (isur.NE.-1) hescal = hescal + flux
            ELSE
              CALL WN('LoadEireneData','Unknown particle type')
            ENDIF          
          ENDDO

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
c slmod begin
        IF (idring(ir).EQ.BOUNDARY) CYCLE
c slmod end        
        DO ik = 1, nks(ir)

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
          pinatom(ik,ir) = pinatom(ik,ir) / kvols(ik,ir)
          pinmol (ik,ir) = pinmol (ik,ir) / kvols(ik,ir)
          pinline(ik,ir,1:6,H_BALPHA) = pinline(ik,ir,1:6,H_BALPHA) /
     .                                  kvols(ik,ir)
          pinline(ik,ir,1:6,H_BGAMMA) = pinline(ik,ir,1:6,H_BGAMMA) /
     .                                  kvols(ik,ir)

c          DO i1 = 1, 6
c            pinline(ik,ir,i1,H_BALPHA) = pinline(ik,ir,i1,H_BALPHA) /
c    .                                    kvols(ik,ir)
c            pinline(ik,ir,i1,H_BGAMMA) = pinline(ik,ir,i1,H_BGAMMA) /
c    .                                    kvols(ik,ir)
c          ENDDO
c...
          pinalpha(ik,ir) = pinline(ik,ir,6,H_BALPHA)

          sumion = sumion + pinion(ik,ir) * kvols(ik,ir)
        ENDDO
      ENDDO


      ir = irsep
      DO ik = 1, 10
        WRITE(0,*) 'CHECK:',
     .    pinline(ik,ir,6,H_BALPHA),
     .    pinline(ik,ir,6,H_BGAMMA),
     .    pinline(ik,ir,6,H_BGAMMA) /
     .    pinline(ik,ir,6,H_BALPHA)
      ENDDO


c      WRITE(0,*) 'SUMION:',sumion

      CALL OutputEIRENE(67,'WORKING WITH EIRENE04')

c...  Make sure that the whole EIRENE data file was there:
      IF (.NOT.goodeof)  CALL ER('LoadEireneData','Problem with '//
     .                           'eirene.transfer file',*99)

      DEALLOCATE(tdata)

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
c subroutine: WriteEireneDatFileInfo
c
c
      SUBROUTINE WriteEireneDatFileInfo
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      REAL    GetFlux

      INTEGER fp,ik,ir
      LOGICAL first
      REAL    balance,sumion,sumrec,sumpuf,sumflx,lossc,lossp

      DATA first /.TRUE./
      SAVE

      fp = 98

      OPEN(UNIT=fp,FILE='eirene04.txt',ACCESS='SEQUENTIAL',
     .     STATUS='UNKNOWN',POSITION='APPEND',ERR=98)

c...  need to calculate the total flux to targets, recom, ionisation, balance, etc. ...


      sumion = 0.0  ! *** Also need sources/sinks that are outside the magnetic grid
      sumrec = 0.0
      sumflx = 0.0
      sumpuf = 0.0
      CALL VolInteg(pinion,3,1,nrs,sumion)
      CALL VolInteg(pinrec,3,1,nrs,sumrec)
      DO ir = irsep, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        sumflx = sumflx - GetFlux(IKLO,ir) + GetFlux(IKHI,ir)  
      ENDDO
c      sumion = sumion / eirtorfrac
c      sumrec = sumrec / eirtorfrac
c      sumflx = sumflx / eirtorfrac

      balance = (sumflx + sumrec + sumpuf) / (sumion + hescpd)

      lossc = (hescpd - hescal) / (sumflx + sumrec + sumpuf) * 100.0
      lossp = (         hescal) / (sumflx + sumrec + sumpuf) * 100.0

      IF (first) THEN 
        CALL HD(fp,'  PARTICLE BALANCE PER ITERATION','PARBAL-HD',5,85)
        WRITE(fp,*)
        WRITE(fp,'(4X,3A5,A7,4A10,2A6)')
     .    'ITER',' STEP','SUBI','BAL.','TARFLX','VOLREC','PUFF',
     .    'IONIZ','CORE','PUMP'
        WRITE(fp,'(4X,3A5,A7,4A10,2A6)')
     .    '     ','     ','     ','  ','(s-1)','(s-1)','(s-1)',
     .    '(s-1)','(%)','(%)'
        first = .FALSE.
      ENDIF

      WRITE(fp,'(4X,3I5,F7.3,1P,4E10.2,0P,2F6.1)')
     .  rel_count,rel_step,rel_iter,
     .  balance,sumflx,sumrec,sumpuf,sumion,lossc,lossp

      CLOSE (fp)

c...  Contribution to .rel.raw (awkward system):
      WRITE(79,'(A,4I6,F6.3,1P,4E10.2,0P,2F6.1)')
     .  '''PARBALAN  1.00''',rel_step,rel_iter,rel_count,7,    ! '7' indicates the number of quantites
     .  balance,sumflx,sumrec,sumpuf,sumion,lossc,lossp

      RETURN
 98   CALL ER('WriteEireneDatFileInfo','Unable to open .txt file',*99)
 99   STOP
      END
c
c ======================================================================
c
c subroutine: ReadParticleTracks
c
c
      SUBROUTINE ReadParticleTracks_04
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER fp,i1,idum1,idum2,count
      REAL    rdum1,rdum2,rdum3

      WRITE(PINOUT,*)
      WRITE(PINOUT,'(A)') ' Reading EIRENE particle tracks:'

      hwalks(1,1) = 999.0
      hwalks(1,2) = 999.0

      fp = 98

      OPEN(UNIT=fp,FILE='eirtrac.dat',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)

      count = 0
      i1    = 0
      DO WHILE (.TRUE..AND.i1.LT.MAXNWS-2)
        READ(fp,*,END=10,ERR=97) idum1,rdum1,rdum2,rdum3,idum2

        IF (idum2.EQ.0.AND.i1.GT.0) THEN
          count = count + 1

          i1 = i1 + 1
          hwalks(i1,1) = HI
          hwalks(i1,2) = HI
        ENDIF
          
        i1 = i1 + 1
        hwalks(i1,1) = rdum1 * 0.01
        hwalks(i1,2) = rdum2 * 0.01 
      ENDDO
10    CONTINUE

      i1 = i1 + 1
      hwalks(i1,1) = 999.0
      hwalks(i1,2) = 999.0

      CLOSE(fp)

      WRITE(PINOUT,'(A,I8)') '   NO. OF TRACKS= ',count
c      WRITE(0     ,'(A,I8)') '   NO. OF TRACKS= ',count      

      RETURN
96    CALL ER('ReadParticleTracks','EOF'     ,*99)
97    CALL ER('ReadParticleTracks','Problems',*99)
c...  The data file could be missing because either the EIRENE run script
c     needs to be updated, or there were not enough particle tracks followed
c     in EIRENE to reach the track number range set for output in the 
c     EIRENE input file (this is set in block 11 in the template file at the
c     moment):
98    CALL WN('ReadParticleTracks','Unable to open data file')
      RETURN
99    STOP
      END


c
c ======================================================================
c ======================================================================
c ======================================================================
c

c subroutine: DumpTriangles
c
      SUBROUTINE DumpTriangles
      USE mod_eirene04
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'


      REAL       TOL
      PARAMETER (TOL=1.0E-07)

      REAL GetMach,GetJsat,GetFlux 

      INTEGER fp,fp2,ik1,ik2,ir1,ir2,i1,i2,trin,id,nl,v1,v1n,v2,v2n,
     .        vern,in,tarside,region
      LOGICAL status
      REAL    x1,x2,y1,y2,deltax,deltay,Bfrac,fact,
     .        tarte,tarti,tarne,tarv,tarflux,tarisat,tarM
      REAL*8  Bx,By,Bz,beta,brat

      INTEGER, ALLOCATABLE :: nlink1(:,:)    ,nlink2(:,:),
     .                        llink1(:,:,:,:),llink2(:,:,:,:),
     .                        triik(:),triir(:),trisf(:),
     .                        trimap(:,:),triside(:,:),trivert(:,:),
     .                        trisurface(:,:)
      REAL, ALLOCATABLE :: trix(:,:),verx(:),verz(:),
     .                     triy(:,:),very(:)


c      RETURN

      WRITE(0,*) 'DUMPING TRIANGLES'

      ALLOCATE(nlink1(MAXNKS,MAXNRS))
      ALLOCATE(nlink2(MAXNKS,MAXNRS))
      ALLOCATE(llink1(MAXNKS,MAXNRS,100,2))
      ALLOCATE(llink2(MAXNKS,MAXNRS,100,2))

      ALLOCATE(triik(MAXNKS*MAXNRS))
      ALLOCATE(triir(MAXNKS*MAXNRS))
      ALLOCATE(trisf(MAXNKS*MAXNRS))
      ALLOCATE(trix(MAXNKS*MAXNRS,3))
      ALLOCATE(triy(MAXNKS*MAXNRS,3))
      ALLOCATE(trimap(MAXNKS*MAXNRS,3))
      ALLOCATE(triside(MAXNKS*MAXNRS,3))
      ALLOCATE(trisurface(MAXNKS*MAXNRS,3))
      ALLOCATE(trivert(MAXNKS*MAXNRS,3))

      ALLOCATE(verx(MAXNKS*MAXNRS))
      ALLOCATE(very(MAXNKS*MAXNRS))
      ALLOCATE(verz(MAXNKS*MAXNRS))

      vern = 0
      verx = 0.0
      very = 0.0
      verz = 0.0


      fp = 99
      fp2 = 98


c...  Build a list of reference quadrilaterals, and a list of connection
c     map references to each cell:    

      nlink1 = 0
      nlink2 = 0
      DO ir1 = 2, nrs
        IF (idring(ir1).EQ.-1) CYCLE
        DO ik1 = 1, nks(ir1)
          DO ir2 = 1, nrs
            IF (idring(ir2).EQ.BOUNDARY) CYCLE
            DO ik2 = 1, nks(ir2)
              IF (ir2.LT.irsep.AND.ik2.EQ.nks(ir2)) CYCLE
              IF     (irouts(ik2,ir2).EQ.ir1.AND.
     .                ikouts(ik2,ir2).EQ.ik1) THEN

                nlink1(ik1,ir1) = nlink1(ik1,ir1) + 1
                llink1(ik1,ir1,nlink1(ik1,ir1),1) = ik2
                llink1(ik1,ir1,nlink1(ik1,ir1),2) = ir2
              ELSEIF (irins (ik2,ir2).EQ.ir1.AND.
     .                ikins (ik2,ir2).EQ.ik1) THEN
                nlink2(ik1,ir1) = nlink2(ik1,ir1) + 1
                llink2(ik1,ir1,nlink2(ik1,ir1),1) = ik2
                llink2(ik1,ir1,nlink2(ik1,ir1),2) = ir2
              ENDIF
            ENDDO
          ENDDO
c...      
          IF (nlink1(ik1,ir1).LE.1) nlink1(ik1,ir1) = 0
          IF (nlink2(ik1,ir1).LE.1) nlink2(ik1,ir1) = 0
c...      Identify cells that adjacent to boundary rings:
          IF (idring(irins (ik1,ir1)).EQ.-1) nlink1(ik1,ir1) = -1
          IF (idring(irouts(ik1,ir1)).EQ.-1) nlink2(ik1,ir1) = -1
        ENDDO
c        DO ik1 = 1, nks(ir1)     
c         IF (ir1.EQ.irsep.AND.ik1.EQ.25) THEN
c           WRITE(0,*)  'DUMP:',ik1,nlink1(ik1,ir1),nlink2(ik1,ir1)
c         ENDIF
c        ENDDO
      ENDDO 

c...  Construct 2 triangles from each regular OEDGE cell:
c
c     counter-clockwise
c     Triangle 1 contains OEDGE cell verticies 1 and 2, and 
c     triangle 2 contains verticies 3 and 4.
c
c
      DO ir1 = 2, nrs
        IF (idring(ir1).EQ.-1) CYCLE

c        IF (ir1.NE.irtrap+1.AND.ir1.NE.irtrap+2) CYCLE
c        IF (ir1.NE.irsep.AND.ir1.NE.nrs.AND.ir1.NE.irsep+1.AND.
c     .      ir1.NE.irsep-1) CYCLE
c        IF (ir1.LT.irtrap) CYCLE

        DO ik1 = 1, nks(ir1)
          IF (ir1.LT.irsep.AND.ik1.EQ.nks(ir1)) CYCLE
 
          id = korpg(ik1,ir1)

          trin = trin + 1
          triik(trin) = ik1
          triir(trin) = ir1
          trisf(trin) = 23
          trix(trin,1) = rvertp(1,id)         
          triy(trin,1) = zvertp(1,id) 
          trix(trin,2) = rvertp(3,id) 
          triy(trin,2) = zvertp(3,id) 
          trix(trin,3) = rvertp(2,id) 
          triy(trin,3) = zvertp(2,id) 

          trin = trin + 1
          triik(trin) = ik1
          triir(trin) = ir1
          trisf(trin) = 14
          trix(trin,1) = rvertp(3,id) 
          triy(trin,1) = zvertp(3,id) 
          trix(trin,2) = rvertp(1,id) 
          triy(trin,2) = zvertp(1,id) 
          trix(trin,3) = rvertp(4,id) 
          triy(trin,3) = zvertp(4,id) 
        ENDDO
      ENDDO

c.... Slice triangles with more than one neighbour:
      i1 = 0
      status = .TRUE.
      DO WHILE (status)
        i1 = i1 + 1

        ik1 = triik(i1)
        ir1 = triir(i1)

        IF     (trisf(i1).EQ.23) THEN
          IF (nlink2(ik1,ir1).LE.0) CYCLE

          nl = nlink2(ik1,ir1)

c...      Make room for new triangles:
          DO i2 = trin, i1+1, -1
            triik(i2+nl-1)   = triik(i2)
            triir(i2+nl-1)   = triir(i2)
            trisf(i2+nl-1)   = trisf(i2)
            trix (i2+nl-1,1) = trix (i2,1)
            triy (i2+nl-1,1) = triy (i2,1)
            trix (i2+nl-1,2) = trix (i2,2)
            triy (i2+nl-1,2) = triy (i2,2)
            trix (i2+nl-1,3) = trix (i2,3)
            triy (i2+nl-1,3) = triy (i2,3)
          ENDDO
          trin = trin + nl - 1
  
c...      Duplicate triangle:
          DO i2 = i1+1, i1+nl-1
            triik(i2) =  triik(i1)
            triir(i2) =  triir(i1)
            trisf(i2) = -trisf(i1)
            trix(i2,1) = trix(i1,1)
            triy(i2,1) = triy(i1,1)
            trix(i2,2) = trix(i1,2)
            triy(i2,2) = triy(i1,2)
            trix(i2,3) = trix(i1,3)
            triy(i2,3) = triy(i1,3)
          ENDDO

c...      Update trailing vertex from connection map data:
          DO i2 = 2, nl
            ik2 = llink2(ik1,ir1,i2,1)
            ir2 = llink2(ik1,ir1,i2,2)
            id = korpg(ik2,ir2)
            trix(i1+i2-1,3) = rvertp(1,id) 
            triy(i1+i2-1,3) = zvertp(1,id) 
          ENDDO

c...      Make triangles in this series consistent with each subsequent
c         triangle in the series:
          DO i2 = 1, nl-1
            trix(i1+i2-1,2) = trix(i1+i2,3)  
            triy(i1+i2-1,2) = triy(i1+i2,3) 
          ENDDO

        ELSEIF (trisf(i1).EQ.14) THEN
          IF (nlink1(ik1,ir1).LE.0) CYCLE

          nl = nlink1(ik1,ir1)

c...      Make room for new triangles:
          DO i2 = trin, i1+1, -1
            triik(i2+nl-1)   = triik(i2)
            triir(i2+nl-1)   = triir(i2)
            trisf(i2+nl-1)   = trisf(i2)
            trix (i2+nl-1,1) = trix (i2,1)
            triy (i2+nl-1,1) = triy (i2,1)
            trix (i2+nl-1,2) = trix (i2,2)
            triy (i2+nl-1,2) = triy (i2,2)
            trix (i2+nl-1,3) = trix (i2,3)
            triy (i2+nl-1,3) = triy (i2,3)
          ENDDO
          trin = trin + nl - 1
  
c...      Duplicate triangle:
          DO i2 = i1+1, i1+nl-1
            triik(i2) =  triik(i1)
            triir(i2) =  triir(i1)
            trisf(i2) = -trisf(i1)
            trix(i2,1) = trix(i1,1)
            triy(i2,1) = triy(i1,1)
            trix(i2,2) = trix(i1,2)
            triy(i2,2) = triy(i1,2)
            trix(i2,3) = trix(i1,3)
            triy(i2,3) = triy(i1,3)
          ENDDO

c...      Update trailing vertex from connection map data:
          DO i2 = 2, nl
            ik2 = llink1(ik1,ir1,i2,1)
            ir2 = llink1(ik1,ir1,i2,2)
            id = korpg(ik2,ir2)
            trix(i1+i2-1,2) = rvertp(2,id) 
            triy(i1+i2-1,2) = zvertp(2,id) 
            IF (ik1.EQ.25.AND.
     .        ir1.EQ.irsep) 
     .        WRITE(0,*) 'ik,ir2=',ik2,ir2,rvertp(2,id),zvertp(2,id),nl
          ENDDO

c...      Make triangles in this series consistent with each subsequent
c         triangle in the series:
          DO i2 = 1, nl-1
            trix(i1+i2-1,3) = trix(i1+i2,2)  
            triy(i1+i2-1,3) = triy(i1+i2,2) 
          ENDDO

        ENDIF

        status = .FALSE.
        IF (i1.LT.trin) status = .TRUE.
      ENDDO


c         DO i2 = 1, trin
c           WRITE(0,'(A,I4,3F10.5)') '->',i2,trix(i2,1),triy(i2,1)
c           WRITE(0,'(A,I4,3F10.5)') '  ',i2,trix(i2,2),triy(i2,2)
c           WRITE(0,'(A,I4,3F10.5)') '  ',i2,trix(i2,3),triy(i2,3)
c         ENDDO
 


c...  Build connection map:

c...  THIS CAN BE MORE EFFICIENT:
      DO i1 = 1, trin
        DO v1 = 1, 3
          v1n = v1 + 1
          IF (v1n.EQ.4) v1n = 1
          DO i2 = 1, trin
            IF (i1.EQ.i2) CYCLE
            DO v2 = 1, 3    
              v2n = v2 - 1
              IF (v2n.EQ.0) v2n = 3
             
              IF (ABS(trix(i1,v1 )-trix(i2,v2 )).LT.TOL.AND.
     .            ABS(triy(i1,v1 )-triy(i2,v2 )).LT.TOL.AND.
     .            ABS(trix(i1,v1n)-trix(i2,v2n)).LT.TOL.AND.
     .            ABS(triy(i1,v1n)-triy(i2,v2n)).LT.TOL) THEN
                trimap (i1,v1) = i2
                triside(i1,v1) = v2n
              ENDIF

            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO i2 = 1, trin
        DO v1 = 1, 3
          IF (trimap(i2,v1).EQ.0) THEN
c...        Identify neighbourless surfaces and assign the appropriate mapping
c           to the appropriate non-default standard surface:
            ik1 = triik(i2)
            ir1 = triir(i2)
            IF     (v1.EQ.3.AND.ik1.EQ.1       .AND.ir1.GE.irsep) THEN
c...          Low IK index target segment, EIRENE poloidal suface 4:
              trisurface(i2,v1) = 2                            
c              trisurface(i2,v1) = 4                            
            ELSEIF (v1.EQ.3.AND.ik1.EQ.nks(ir1).AND.ir1.GE.irsep) THEN
c...          High IK index target segment, EIRENE poloidal surface 5:
              trisurface(i2,v1) = 3
c              trisurface(i2,v1) = 5                            
            ELSEIF (v1.EQ.2.AND.irins(ik1,ir1) .EQ.1) THEN
c...          Core boundary, EIRENE radial surface 1:
              trisurface(i2,v1) = 1                            
            ELSEIF ((v1.EQ.2.AND.irouts(ik1,ir1).EQ.irwall).OR.
c...          Added for MAST-DN secondary PFZ:
     .              (v1.EQ.2.AND.irins (ik1,ir1).EQ.irwall)) THEN
c...          SOL boundary, EIRENE radial surface 4:
              trisurface(i2,v1) = 4                            
c              trisurface(i2,v1) = 2                            
            ELSEIF (v1.EQ.2.AND.irins(ik1,ir1) .EQ.irtrap) THEN
c...          PFZ boundary, EIRENE radial surface 5:
              trisurface(i2,v1) = 5
c              trisurface(i2,v1) = 3                            
            ELSEIF (v1.EQ.1.OR.v1.EQ.3) THEN
              WRITE(0,*) 'PROBLEM 3:',i2,v1,triik(i2),triir(i2)
            ENDIF
             
          ELSEIF (v1.EQ.2.AND.ABS(trisf(i2)).EQ.23.AND.
     .            (triir(i2).EQ.irsep-1.OR.triir(i2).EQ.nrs)) THEN
c...        Label "right" side of separatrix as EIRENE non-standard surface 6:
            trisurface(i2,v1) = 6
          ELSEIF (v1.EQ.2.AND.ABS(trisf(i2)).EQ.14.AND.
     .            triir(i2).EQ.irsep) THEN
c...        Label "left" side of separatrix as EIRENE non-standard surface 7:
            trisurface(i2,v1) = 7

          ELSEIF (trimap(trimap(i2,v1),triside(i2,v1)).NE.i2) THEN
            WRITE(0,*) 'PROBLEM 4:',i2,v1,triik(i2),triir(i2)

          ELSEIF (triside(trimap(i2,v1),triside(i2,v1)).NE.v1) THEN
            WRITE(0,*) 'PROBLEM 5:',i2,v1,triik(i2),triir(i2)

          ENDIF
        ENDDO
      ENDDO

c...  All is well, so generate a master list of points with pointers to the respective
c     triangle verticies:
      DO i1 = 1, trin
        DO v1 = 1, 3

c...      Scan list of verticies and assign the relevant vertex
c         number to the triangle:
          status = .FALSE.
          DO i2 = 1, vern
            IF (ABS(verx(i2)-trix(i1,v1)).LT.TOL.AND.
     .          ABS(very(i2)-triy(i1,v1)).LT.TOL) THEN
              IF (trivert(i1,v1).EQ.0) THEN
                trivert(i1,v1) = i2
                status = .TRUE.
              ELSE
                CALL ER('DumpTriangles','Redundant vertex',*99)
              ENDIF
            ENDIF
          ENDDO
c...      Vertex was not found, so add it to the list:
          IF (.NOT.status) THEN
            vern = vern + 1
            verx(vern) = trix(i1,v1)
            very(vern) = triy(i1,v1)
            trivert(i1,v1) = vern
          ENDIF

        ENDDO
      ENDDO

c      WRITE(0,*) 'vern:',vern,trin

c...  Dump triangles:
      OPEN(UNIT=fp,FILE='triangles.dat',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) trin
      DO i1 = 1, trin
        WRITE(fp,'(I6,3(2F10.6,2X))')
     .    i1,(verx(trivert(i1,i2)),very(trivert(i1,i2)),i2=1,3)
c        WRITE(fp,'(I6,3(2F10.6,2X))')
c     .    i1,(trix(i1,i2),triy(i1,i2),i2=1,3)
      ENDDO
      CLOSE(fp)      

c...  Dump triangles:
      OPEN(UNIT=fp,FILE='triangles.points',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) vern
      DO i1 = 1, vern
        WRITE(fp,'(I6,3F12.6)') i1,verx(i1)*100.0,very(i1)*100.0,0.0
      ENDDO
      CLOSE(fp)      

c...  Dump triangles:
      OPEN(UNIT=fp,FILE='triangles.sides',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) trin
      DO i1 = 1, trin
        WRITE(fp,'(I6,4X,3I6)') i1,(trivert(i1,v1),v1=1,3)
      ENDDO
      CLOSE(fp)      

c...  Dump triangles:
      OPEN(UNIT=fp,FILE='triangles.map',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) trin
      DO i1 = 1, trin
        WRITE(fp,'(I6,4X,3(3I6,4X),2I6)') i1,
     .    (trimap(i1,v1),triside(i1,v1),trisurface(i1,v1),v1=1,3),
     .    triik(i1),triir(i1)
      ENDDO
      CLOSE(fp)      


      fact = qtim * qtim * emi / crmi

c...  Dump plasma data:
      OPEN(UNIT=fp ,FILE='triangles.plasma',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      OPEN(UNIT=fp2,FILE='triangles.efield',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      

c...  Header:
      WRITE(fp,'(A)') '* VERSION 1.0 OSM PLASMA FILE FOR TRIANGULAR '//
     .                'GRID'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* BULK PLASMA DATA'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A7,2A8,2X,A10,2X,3A10,2X,3A12)') 
     .  '* Index','Te','Ti','ne','vx',
     .  'vy','vz','Bx','By','Bz'
      WRITE(fp,'(A7,2A8,2X,A10,2X,3A10,2X,3A12)') 
     .  '*      ','(eV)','(eV)','(cm-3)','(cm s-1)',
     .  '(cm s-1)','(cm s-1)','(Tesla)','(Tesla)','(Tesla)'


c...  Header:
      WRITE(fp2,'(A)') '* VERSION 1.0 OSM PLASMA FILE FOR TRIANGULAR '//
     .                'GRID'
      WRITE(fp2,'(A)') '*'
      WRITE(fp2,'(A)') '* BULK PLASMA DATA (PART DEUX)'
      WRITE(fp2,'(A)') '*'
      WRITE(fp2,'(A7,3A12)') 
     .  '* Index','Ex','Ey','Ez'
      WRITE(fp2,'(A7,3A12)')
     .  '*      ','V m-1','V m-1','V m-1'

      
      WRITE(fp,*) trin
      DO i1 = 1, trin
        ik1 = triik(i1)
        ir1 = triir(i1)
        id = korpg(ik1,ir1)
        x1 = 0.5 * (rvertp(1,id) + rvertp(2,id))
        y1 = 0.5 * (zvertp(1,id) + zvertp(2,id))
        x2 = 0.5 * (rvertp(3,id) + rvertp(4,id))
        y2 = 0.5 * (zvertp(3,id) + zvertp(4,id))
        deltax = (x2 - x1)
        deltay = (y2 - y1)
        brat = DBLE(bratio(ik1,ir1))
        beta = DBLE(deltax / deltay)
        Bz = DSQRT(1.0 - brat**2)
        By = brat * DSQRT(1.0D0 / (1.0D0 + beta**2)) * 
     .       DBLE(SIGN(1.0,deltay))
        Bx = beta * By
c...    CBPHI is the on-axis B-field value specified in the OSM input file:
        Bfrac = cbphi * r0 / rs(ik1,ir1) 

        WRITE(fp,'(I7,2F8.2,2X,1P,E10.2,2X,3E10.2,2X,
     .             3E12.4,0P,6X,2I4)') i1,
c     .             3E12.4,0P,4X,2I6,3F10.5)') i1,
c        WRITE(fp,'(I6,2F10.2,1P,E10.2,2X,5E10.2,0P,2I6,2F10.2))') i1,
     .         ktebs(ik1,ir1),ktibs(ik1,ir1),knbs(ik1,ir1)*1.0E-06,
     .         SNGL(Bx)*kvhs(ik1,ir1)*100.0,  !/qtim,
     .         SNGL(By)*kvhs(ik1,ir1)*100.0,  !/qtim,
     .         SNGL(Bz)*kvhs(ik1,ir1)*100.0,  !/qtim,
     .         SNGL(Bx)*Bfrac,SNGL(By)*Bfrac,SNGL(Bz)*Bfrac,
c...    Target quantities:
c     .         tarside,tarflux,tarte,tarti,tarne*1.0E-06,tarv*100.0,
c     .         tarM,tarisat*1.0E-04, 
c...    Grid cell where triangle originated
     .         ik1,ir1
c     .         SQRT(SNGL(Bx)**2+
c     .              SNGL(By)**2+
c     .              SNGL(Bz)**2),bratio(ik1,ir1),
c     .         SQRT(SNGL(Bx)**2+
c     .              SNGL(By)**2)
c     .         SQRT((SNGL(Bx)*kvhs(ik1,ir1))**2+
c     .              (SNGL(By)*kvhs(ik1,ir1))**2+
c     .              (SNGL(Bz)*kvhs(ik1,ir1))**2)*
c     .              sign(1.0,kvhs(ik1,ir1)),
c     .         kvhs(ik1,ir1),
c     .         ik1,ir1,
c     .         bratio(ik1,ir1),SNGL(DSQRT(Bx**2+By**2)/Bz)
cc     .         SNGL(Bx),SNGL(By),SNGL(Bz),
cc     .         DSQRT(Bx**2+By**2+Bz**2),ik1,ir1,
cc     .         bratio(ik1,ir1),SNGL(DSQRT(Bx**2+By**2)/Bz)


c...    Dump efield data:
        WRITE(fp2,'(I7,1P,3E12.4,10X,0P,2I4,1P,E12.4)') i1,
     .         SNGL(Bx)*kes(ik1,ir1)/fact,
     .         SNGL(By)*kes(ik1,ir1)/fact,
     .         SNGL(Bz)*kes(ik1,ir1)/fact,
     .         ik1,ir1,kes(ik1,ir1)/fact

      ENDDO

c      WRITE(fp,'(A)') '*'
c      WRITE(fp,'(A)') '* TARGET DATA'
      WRITE(fp,'(A)') '*'
c      WRITE(fp,'(A7,A6,A10,2A8,2A10,A6,A10') 
c     .  '* Index','Side','flux','Te','Ti','ne','v','Mach','Isat'
c      WRITE(fp,'(A7,A6,A10,2A8,2A10,A6,A10') 
c     .  '*      ','','(A)','(eV)','(eV)','(cm-3)','(cm s-1)','',
c     .  '(A cm-2)'
      
      WRITE(fp,*) 2*(nrs-irsep-1)

      DO i1 = 1, trin
        ik1 = triik(i1)
        ir1 = triir(i1)
        IF (ik1.NE.1.AND.ik1.NE.nks(ir1).OR.ir1.LT.irsep) CYCLE
        DO v1 = 1, 3
          IF (v1.NE.3.OR.trisurface(i1,v1).EQ.0) CYCLE
c...      Process target data:
          IF (ik1.EQ.1) THEN
            in = idds(ir1,2)
            region = IKLO
          ELSE
            in = idds(ir1,1)
            region = IKHI
          ENDIF
          tarside = v1
          tarte = kteds(in)
          tarti = ktids(in)
          tarne = knds(in)
          tarv = kvds(in)
          tarM = GetMach(tarv,tarte,tarti)
          tarisat = GetJsat(tarte,tarti,tarne,tarv)
          tarflux = GetFlux(region,ir1) / eirtorfrac / eirsrcmul * ECH

          WRITE(fp,'(I7,I6,1P,E10.2,0P,2F8.2,1P,2E10.2,0P,
     .               F6.2,1P,E10.2,0P,6X,2I4)') 
c...    Target quantities:
     .           i1,tarside,tarflux,tarte,tarti,tarne*1.0E-06,
     .           tarv*100.0,tarM,tarisat*1.0E-04, 
c...    Grid cell where triangle originated
     .           ik1,ir1
        ENDDO
      ENDDO

c...  All done:
      CLOSE(fp)      
      CLOSE(fp2)      

c...  Store triangle list in case opacity data is read from a separate
c     EIRENE run:
      IF (trin.GT.MAXNKS*MAXNRS) 
     .  CALL ER('DumpTrinagles','Arrays bounds error',*99)
      eirtrin = trin
      DO i1 = 1, trin
        eirtriik(i1) = triik(i1)
        eirtriir(i1) = triir(i1)
      ENDDO

c...  mod_eirene04:
      ntri = trin
      nver = vern

      DEALLOCATE(nlink1)
      DEALLOCATE(nlink2)
      DEALLOCATE(llink1)
      DEALLOCATE(llink2)

      DEALLOCATE(triik)
      DEALLOCATE(triir)
      DEALLOCATE(trix)
      DEALLOCATE(triy)
      DEALLOCATE(trimap)
      DEALLOCATE(triside)
      DEALLOCATE(trisurface)
      DEALLOCATE(trivert)

      DEALLOCATE(verx)
      DEALLOCATE(very)
      DEALLOCATE(verz)

c      STOP 'DUMPING'
      WRITE(0,*) 'DONE'

      RETURN
96    WRITE(0,*) 'DUMPTRIANGES: PROBLEMS WITH FILE ACCESS'
      STOP
99    STOP
      END
