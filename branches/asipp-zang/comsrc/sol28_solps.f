c     -*-Fortran-*-
c
c ======================================================================
c
c subroutine: SetupToroidalSurfaces
c
c
      SUBROUTINE ReadSOLPSDataFile_OSM(idata)
      USE mod_solps_params
      USE mod_solps
      IMPLICIT none

      INTEGER, INTENT(IN) :: idata
 
      INTEGER   fp,outfp,column,idum(2),maxik,maxir,count
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

c...  First time the routine is called -- as indicated by the 
c     unallocated MAP_DIVIMP array -- so sniff the file and 
c     see how large the SOLPS data arrays need to be:
      IF (.NOT.ALLOCATED(solps_ik)) THEN
        first_call = .TRUE.
        maxik = 0
        maxir = 0
        SELECTCASE (solps_data(idata)%format)
          CASE (1)  ! AK format for i1541
            READ(fp,*)
            READ(fp,*)
            READ(fp,*)
            DO WHILE(.TRUE.)
              READ(fp,*,END=10,ERR=91) idum(1:2)
              maxik = MAX(maxik,idum(1))
              maxir = MAX(maxir,idum(2))
            ENDDO
 10         CONTINUE
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
        CASE (1) 
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
              solps_ik (count)     = idum(1)
              solps_ir (count)     = idum(2)
              solps_cen(count,1:2) = rdum(1:2) ! * 1.0D-03 Rescaling removed on 18/08/2010 after noticing that AK modified 
            ENDIF                              !           the 2dvi script for extracting SOLPS data. -SL
          ENDDO
 50       CONTINUE
        CASE DEFAULT
          CALL ER('ReadSOLPSDataFile_OSM','Unknown format',*99)
      ENDSELECT

      CLOSE(fp)



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
      SUBROUTINE LoadSOLPSData_OSM
      USE mod_sol28_global
      USE mod_solps
      USE mod_grid_divimp
      IMPLICIT none

      CHARACTER*16, PARAMETER ::
     .  type_name(4) = ['density','velocity','pressure','temperaure']

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
            DO i = 1, solps_maxik*solps_maxir
c...          SONNETIK,IR, which map the DIVIMP grid cells to the SONNET
c             grid indeces, were assigned when the grid was read in: 
              IF (solps_ik(i).EQ.divimp_ik(ik,ir).AND.  
     .            solps_ir(i).EQ.divimp_ir(ik,ir)) THEN
                IF (map_divimp(ik,ir).NE.0)
     .            CALL ER('LoadSOLPSData','MAP_DIVIMP error',*99)
                map_divimp(ik,ir) = i
                EXIT
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF


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
