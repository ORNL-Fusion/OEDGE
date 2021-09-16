c     -*-Fortran-*-
c
c ======================================================================
c
c subroutine: ReadSOLPSDataFile_OSM(idata)
c
c
      SUBROUTINE ReadSOLPSDataFile_OSM(idata)
      USE mod_solps_params
      USE mod_solps
      IMPLICIT none

      INTEGER, INTENT(IN) :: idata
 
      INTEGER   fp,outfp,column,idum(2),maxik,maxir,count,i,shift
      LOGICAL   output,first_call
      CHARACTER filename*1024,buffer*1024
      REAL      rdum(10)

      SAVE

      first_call = .FALSE.
      output = .TRUE.
      outfp = 0
      fp = 99

      filename = TRIM(solps_data(idata)%fname)

      WRITE(outfp,*) 'LOADING SOLPS DATA '//TRIM(filename)

      OPEN(UNIT=fp,FILE=TRIM(filename),
     .     ACCESS='SEQUENTIAL',STATUS='OLD',ERR=90)

      shift = 0
      IF (solps_indexing.EQ.1) shift = -1
      IF (solps_indexing.EQ.2) shift =  1

c...  First time the routine is called -- as indicated by the 
c     unallocated MAP_DIVIMP array -- so sniff the file and 
c     see how large the SOLPS data arrays need to be:
      IF (.NOT.ALLOCATED(solps_ik)) THEN
        first_call = .TRUE.
        maxik = 0
        maxir = 0
        count = 0
        SELECTCASE (solps_data(idata)%format)
          CASE (1:2)  ! AK format for i1541
            READ(fp,*)
            READ(fp,*)
            READ(fp,*)
            DO WHILE(.TRUE.)
              READ(fp,*,END=10,ERR=91) idum(1:2)
              maxik = MAX(maxik,idum(1))
              maxir = MAX(maxir,idum(2))
              count = count + 1
            ENDDO
 10         CONTINUE
            IF (maxik.EQ.1.AND.maxir.EQ.1) maxik = count  ! For Yannick's crap data files
            maxik = maxik + shift
            maxir = maxir + shift
            ALLOCATE(solps_ik (maxik*maxir)) 
            ALLOCATE(solps_ir (maxik*maxir)) 
            ALLOCATE(solps_cen(maxik*maxir,2)) 
          CASE DEFAULT
            CALL ER('ReadSOLPSDataFile_OSM','Unknown format',*99)
        ENDSELECT
        IF (output) WRITE(outfp,*) 'MAXIK,MAXIR=',maxik,maxir
        solps_maxik = maxik
        solps_maxir = maxir
        REWIND(fp)  ! Reset the file pointer
      ENDIF

c...  Load the data:
      maxik = solps_maxik
      maxir = solps_maxir

      ALLOCATE(solps_data(idata)%data(maxik*maxir))
      solps_data(idata)%data = 0.0

      column = solps_data(idata)%column
      count = 0
 
      SELECTCASE (solps_data(idata)%format)
        CASE (1:2) 
          READ(fp,*)
          READ(fp,*)
          READ(fp,*)
          DO WHILE(.TRUE.)
            READ(fp,*,END=50,ERR=91) idum(1:2),rdum(1:column+2)
            IF (idum(1).GT.maxik.OR.idum(2).GT.maxir)
     .        CALL ER('ReadSOLPSDataFile_OSM','Index inconsistency '// 
     .                'discovered when reading SOLPS data',*99)       
            count = count + 1
            solps_data(idata)%data(count) = rdum(column+2)
            IF (first_call) THEN
              solps_n              = count
              solps_ik (count)     = idum(1) + shift
              solps_ir (count)     = idum(2) + shift
              solps_cen(count,1:2) = rdum(1:2) ! * 1.0D-03 Rescaling removed on 18/08/2010 after noticing that AK modified 
            ENDIF                              !           the 2dvi script for extracting SOLPS data. -SL
          ENDDO
 50       CONTINUE
        CASE DEFAULT
          CALL ER('ReadSOLPSDataFile_OSM','Unknown format',*99)
      ENDSELECT

      CLOSE(fp)


c     Rescale the cell positions if they are in mm (older data):
      IF (first_call) THEN
        i = MAXLOC(solps_cen(1:count,1),1)
        IF (solps_cen(i,1).GT.30.0) THEN   ! Hopefully no future machine will be larger than 30 m...
          CALL WN('ReadSOLPSDataFile_OSM','Rescaling SOLPS '//
     .            'coordinate data to metres from millimetres')
          solps_cen = solps_cen * 1.0E-3
        ENDIF
      ENDIF


      RETURN
 90   CALL ER('LoadSOLPSData','File not found',*99)
 91   CALL ER('LoadSOLPSData','Error reading file',*99)
 99   WRITE(0,*) ' FILENAME=',TRIM(filename)
      WRITE(0,*) ' IDATA   =',idata           
      WRITE(0,*) ' MAXIK,IR=',maxik,maxir     
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE solps_AssignDIVIMPIndex(ik,ir,is)
      USE     mod_solps
      use mod_params
      use mod_cgeom
      use mod_cedge2d
      use mod_comtor
      IMPLICIT none
c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'cedge2d'
c     INCLUDE 'comtor'


      INTEGER, INTENT(IN ) :: ik,ir
      INTEGER, INTENT(OUT) :: is


      INTEGER, PARAMETER :: NUMZONE = 5

      INTEGER i,j,n,iz,ic,izone(NUMZONE*NUMZONE+1),outfp
      LOGICAL first_call,debug
      REAL    vrmin,vzmin,vrmax,vzmax,rspan,zspan,r,z,dist,dist_min

      INTEGER, ALLOCATABLE :: zone(:),lzone(:)

      DATA first_call /.TRUE./
      SAVE

      debug = .TRUE.
      outfp = 88

c...  Clean arrays:
      IF (ik.EQ.-1) THEN 
        IF (ALLOCATED( zone)) DEALLOCATE( zone)
        IF (ALLOCATED(lzone)) DEALLOCATE(lzone)
        RETURN
      ENDIF

      n = solps_maxik * solps_maxir

      IF (first_call) THEN
        first_call = .FALSE.

        ALLOCATE( zone(n))
        ALLOCATE(lzone(n))

c...    Assign cell zone assignment based on its spatial location in the
c       grid.  This (dramatically) increases the speed of the search
c       algorithm when building the connection map (grid structure) below.
        vrmin =  HI
        vrmax = -HI
        vzmin =  HI
        vzmax = -HI
c...    Find the spatial extent of the grid:
        DO i = 1, n      
          IF (solps_cen(i,1).LT.vrmin) vrmin = solps_cen(i,1)
          IF (solps_cen(i,1).GT.vrmax) vrmax = solps_cen(i,1)
          IF (solps_cen(i,2).LT.vzmin) vzmin = solps_cen(i,2)
          IF (solps_cen(i,2).GT.vzmax) vzmax = solps_cen(i,2)
        ENDDO
        vrmin = vrmin - 0.001D0  ! Expand the domain slightly...
        vrmax = vrmax + 0.001D0
        vzmin = vzmin - 0.001D0
        vzmax = vzmax + 0.001D0
c...    Assign each cell to a zone:
        rspan = (vrmax - vrmin) / DBLE(NUMZONE)
        zspan = (vzmax - vzmin) / DBLE(NUMZONE)
        IF (debug) THEN
          WRITE(outfp,*) 'GRID DOMAIN EXTENT:'
          WRITE(outfp,*) '  vrmin',vrmin
          WRITE(outfp,*) '  vrmax',vrmax
          WRITE(outfp,*) '  vzmin',vzmin
          WRITE(outfp,*) '  vzmax',vzmax
          WRITE(outfp,*) '  rspan',rspan
          WRITE(outfp,*) '  zspan',zspan
        ENDIF
        DO i = 1, n
          zone(i) =  INT( (solps_cen(i,1) - vrmin) / rspan ) + 1 +
     .              (INT( (solps_cen(i,2) - vzmin) / zspan )     ) * 
     .              NUMZONE
c          write(88,*) 'check',INT( (solps_cen(i,1) - vrmin) / rspan )+1,
c     .                        INT( (solps_cen(i,2) - vzmin) / zspan )     
        ENDDO
c...    Create list of zone assigments:
        j = 0
        DO iz = 1, NUMZONE*NUMZONE
          izone(iz) = j + 1
          DO i = 1, n
            IF (zone(i).EQ.iz) THEN 
              j = j + 1
              lzone(j) =  i
               zone(i) = -1
            ENDIF
          ENDDO
        ENDDO
        izone(iz) = j + 1

        DO i = 1, n
          IF (zone(i).NE.-1) THEN 
            STOP 'SHIT BOMBS'
          ENDIF
        ENDDO

c...    
        WRITE(0,*) 'ZONES:'
        WRITE(0,*) izone
        WRITE(0,*) 'n', n
        
      ENDIF

c...  Identify which zone the DIVIMP grid cell is in:
      r = rs(ik,ir)
      z = zs(ik,ir)
      iz = MAX(1,MIN(NUMZONE  ,INT((r-vrmin)/rspan)+1)) +
     .     MAX(0,MIN(NUMZONE-1,INT((z-vzmin)/zspan)  )) * 
     .     NUMZONE

      IF (debug) WRITE(outfp,'(A,3I6,2F14.6,2X,4F14.6,2X,2F14.6)') 
     .             'searching',ik,ir,iz,
     .             r,z,
     .             vrmin,vzmin,vrmax,vzmax,
     .             rspan,zspan

      dist_min = HI
      is       = -1
      DO j = izone(iz), izone(iz+1)-1 
        i    = lzone(j)
        dist = SQRT((r - solps_cen(i,1))**2 + (z - solps_cen(i,2))**2)
        IF (dist.LT.dist_min) THEN
          dist_min = dist
          is       = i
          IF (debug) WRITE(outfp,'(A,2I6,2F14.6,2X,F14.6,I6) ')
     .      'dist_min',ik,ir,r,z,dist_min,i
        ENDIF
      ENDDO

      IF (is.EQ.-1) THEN
        STOP 'dcasmdad -- not found!'
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadSOLPSData_OSM
      USE mod_sol28_global
      USE mod_solps
      USE mod_grid_divimp
      IMPLICIT none

      CHARACTER*10, PARAMETER ::
     .  type_name(4) = ['density   ','velocity  ','pressure  ',
     .                  'temperaure']
c     .  type_name(4) = ['density','velocity','pressure','temperaure']  ! gfortran

      INTEGER i,ik,ir
 
      WRITE(logfp,*) 
      WRITE(logfp,*) 'LOADING SOLPS DATA:'      
      WRITE(logfp,'(A7,2X,A9,13X,3(A6,2X),12X,A3,19X,3(A6,2X))') 
     .  'Version','File Name','Format','Column','Type','Tag','Z','A', 
     .  'Charge'

      DO i = 1, nsolps_data

        WRITE(logfp,'(F7.2,2X,A20,2X,3(I4,4X),A12,A20,2X,3(I4,4X))')
     .    solps_data(i)%version,
     .    solps_data(i)%fname,
     .    solps_data(i)%format,
     .    solps_data(i)%column,
     .    solps_data(i)%type,
     .    type_name(solps_data(i)%type),
     .    solps_data(i)%tag,
     .    solps_data(i)%z,
     .    solps_data(i)%a,
     .    solps_data(i)%charge

c...  
        IF (ALLOCATED(solps_data(i)%data)) 
     .    CALL ER('LoadSOLPSData_OSM','DATA array already '//
     .            'allocated for this item',*99)       
 
        CALL ReadSOLPSDataFile_OSM(i)
      ENDDO

c...  Map the DIVIMP array structure to the SOLPS data array structure:
      IF (ALLOCATED(divimp_ik)) THEN
        IF (ALLOCATED(map_divimp)) DEALLOCATE(map_divimp)
        ALLOCATE(map_divimp(divimp_maxnks,divimp_maxnrs)) 
        map_divimp = 0
        DO ir = 1, divimp_maxnrs
          DO ik = 1, divimp_maxnks
            IF (divimp_ik(ik,ir).EQ.0) CYCLE

            IF (solps_maxir.EQ.1) THEN 

              stop 'what?'

              write(88,*) 'processing next cell -- ik,ir',ik,ir

              CALL solps_AssignDIVIMPIndex(ik,ir,map_divimp(ik,ir))

c              write(0,*) 'processing ik,ir',ik,ir,map_divimp(ik,ir)

            ELSE

              DO i = 1, solps_maxik*solps_maxir
c...            SONNETIK,IR, which map the DIVIMP grid cells to the SONNET
c               grid indices, were assigned when the grid was read in: 
                IF (solps_ik(i).EQ.divimp_ik(ik,ir).AND.  
     .              solps_ir(i).EQ.divimp_ir(ik,ir)) THEN
                  IF (map_divimp(ik,ir).NE.0)
     .              CALL ER('LoadSOLPSData','MAP_DIVIMP error',*99)
                  map_divimp(ik,ir) = i
                  EXIT
                ENDIF
              ENDDO

            ENDIF

          ENDDO
        ENDDO
      ENDIF


c...  Clean allocated arrays:
      CALL solps_AssignDIVIMPIndex(-1,-1,-1)


c      DO ir = 1, solps_maxir
c        DO ik = 1, solps_maxik
c          WRITE(0,*) 'DATA:',ik,ir,
c     .      solps_data(nsolps_data)%data(ik+(ir-1)*solps_maxik)
c        ENDDO
c      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
