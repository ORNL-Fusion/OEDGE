c     -*-Fortran-*-
c
c ======================================================================
c
      LOGICAL FUNCTION osmGetLine(fp,buffer,mode)
      USE mod_sol28_io
      USE mod_sol28_global
      IMPLICIT none

      INTEGER  , INTENT(IN)  :: fp,mode
      CHARACTER, INTENT(OUT) :: buffer*(*)

      INTEGER i,n

      osmGetLine = .TRUE. 

      DO WHILE(.TRUE.)
        IF (fp.EQ.-1) THEN
          iiteration = iiteration + 1
          buffer = iteration_buffer(iiteration)
        ELSE
          READ(fp,'(A)',END=10) buffer
        ENDIF

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

c        WRITE(0,*) '  BUFFER:',TRIM(buffer),mode

        IF (LEN_TRIM(buffer).EQ.0) THEN
c...      Comment or blank line, so continue:
        ELSEIF (mode.EQ.ALL_LINES) THEN
          EXIT
        ELSEIF (buffer(1:1).EQ.'{'.OR.buffer(2:2).EQ.'{') THEN
          IF (mode.EQ.WITH_TAG) THEN
c...        Tag line: 
            EXIT
          ELSE
c...        Tag line was not requested so backup the file position:
            osmGetLine = .FALSE.
            IF (fp.EQ.-1) THEN
              iiteration = iiteration - 1
            ELSE
              BACKSPACE(fp)
            ENDIF
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
c      WRITE(0,*) 'BUFFER:',buffer(1:50),osmGetLine
      RETURN

 10   osmGetLine = .FALSE.
      RETURN
      END
c
c ======================================================================
c
      SUBROUTINE PadBufferI(npad,buffer,padding)
      IMPLICIT none

      INTEGER,   INTENT(IN)  :: npad,padding
      CHARACTER, INTENT(OUT) :: buffer*(*)

      INTEGER i,j

      i = LEN_TRIM(buffer)

      WRITE(buffer(i+2:i+256),*) (padding,j=1,npad)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE PadBufferR(npad,buffer,padding)
      IMPLICIT none

      INTEGER,   INTENT(IN)  :: npad
      CHARACTER, INTENT(OUT) :: buffer*(*)
      REAL   ,   INTENT(IN)  :: padding

      INTEGER i,j

      i = LEN_TRIM(buffer)

      WRITE(buffer(i+2:i+256),*) (padding,j=1,npad)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: ProcessIterationBlocks
c
c Moved to sol28_utility.f (in DIVIMP) for compatibility with OUT.
c
c ======================================================================
c
      SUBROUTINE LoadOptions
      USE mod_sol28_params
      USE mod_sol28_io
      USE mod_sol28_global
      USE mod_legacy
      IMPLICIT none

      LOGICAL osmGetLine

      INTEGER   fp,i
      CHARACTER buffer*1024

      LOGICAL :: status = .TRUE. 

c...  Set default values:
      CALL osm_InitializeOptions

c...  Open log file (may be closed below if logfp is set to 0):
c      OPEN(UNIT=logfp,FILE='osm_log.dat',ACCESS='SEQUENTIAL',
c     .     STATUS='REPLACE',ERR=98)     

c...  Open OSM input file that contains option settings:
      fp = 99
      OPEN(UNIT=fp,FILE='osm.input',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)

c...  Scan input file looking for tags which are marked by a '{' in 
c     column 2 on an input line and contained by curly brackets:
      DO WHILE(osmGetLine(fp,buffer,WITH_TAG))
c...    Isolate tag string:
        DO i = 2, LEN_TRIM(buffer)
          IF (buffer(i:i).EQ.'}') EXIT
        ENDDO
        CALL ProcessInputTag(fp,i,buffer,status)
        IF (.NOT.status) EXIT
      ENDDO
      CLOSE(fp)

      CALL ProcessIterationBlocks


      RETURN
 98   CALL ER('LoadOptions','Error reading file',*99)
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE ProcessInputTag(fp,itag,buffer,status)
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_legacy
      IMPLICIT none

      INTEGER, INTENT(IN) :: fp,itag
      LOGICAL   status
      CHARACTER buffer*(*)

      INTEGER i

      WRITE(logfp,*) 'TAG:>'//buffer(3:itag-1)//'<'

      DO i = 3, itag-1
        IF (buffer(i:i).EQ.' ') EXIT
      ENDDO
      i = i - 1

c...  Check for reserved group markers which are indicated by the
c     first two characters in the tag, otherwise treat the tag as
c     miscellaneous:
      SELECTCASE (buffer(3:i))
c      SELECTCASE (buffer(3:4))
        CASE ('C','CO','CON')
          CALL LoadControlOption(fp,buffer,itag)
        CASE ('D','DIV')
          CALL LoadDivimpOption(fp,buffer,itag)
        CASE ('E','EI','EIR')
          CALL LoadEireneOption(fp,buffer,itag)
        CASE ('FI')
          CALL LoadFilamentOption(fp,buffer,itag)
        CASE ('F','GRID','FILE')
          CALL LoadFileOption(fp,buffer,itag)
        CASE ('S')
          CALL LoadSourceOption(fp,buffer,itag)
        CASE ('SOL')
          CALL LoadSolverOption(fp,buffer,itag)
        CASE ('SOLPS')
          CALL LoadSOLPSOption(fp,buffer,itag)
        CASE ('T')
          CALL LoadTransportOption(fp,buffer,itag)
        CASE ('WALL')
          CALL LoadWallOption(fp,buffer,itag)
        CASE DEFAULT
          CALL LoadMiscOption(fp,buffer,itag,status)
      ENDSELECT


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE ReadOptionI(buffer,noption,option)
      IMPLICIT none

      CHARACTER, INTENT(IN)  :: buffer*(*)
      INTEGER  , INTENT(OUT) :: option(*)
      INTEGER  , INTENT(IN)  :: noption
  
      INTEGER        i1
      CHARACTER*1024 cdum1

      CALL PadBufferI(15,buffer,999)

      READ(buffer,*) cdum1,option(1:noption)

      DO i1 = 2, noption
        IF (option(i1).EQ.999) option(i1) = option(1)
      ENDDO

c      WRITE(0,*) 'OPTION:',option(1:noption)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE ReadOptionR(buffer,noption,option)
      IMPLICIT none

      CHARACTER, INTENT(IN)  :: buffer*(*)
      REAL     , INTENT(OUT) :: option(*)
      INTEGER  , INTENT(IN)  :: noption
  
      INTEGER        i1
      CHARACTER*1024 cdum1

      CALL PadBufferR(15,buffer,999.0)

c      WRITE(0,'(a)') 'BUFFER:'//buffer(1:100)

      READ(buffer,*) cdum1,option(1:noption)

      DO i1 = 2, noption
        IF (option(i1).EQ.999.0) option(i1) = option(1)
      ENDDO

c      WRITE(0,*) 'OPTION:',option(1:noption)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadControlOption(fp,buffer,itag)
      USE mod_sol28_io
      USE mod_sol28_global
      IMPLICIT none

      INTEGER  , INTENT(IN) :: fp,itag
      CHARACTER, INTENT(IN) :: buffer*(*)

      LOGICAL osmGetLine

      INTEGER   i
      LOGICAL   finished
      CHARACTER cdum1*512

      SELECTCASE (buffer(3:itag-1))
        CASE('C KIN')
          READ(buffer,*) cdum1,opt%nflukin
        CASE('C EIR')
          READ(buffer,*) cdum1,opt%eirene
        CASE('C EIR_TIME')
          READ(buffer,*) cdum1,opt_eir%time
        CASE('CON EIRENE NTIME')
          CALL ReadOptionI(buffer,1,opt_eir%ntime) 
        CASE('CON EIRENE DTIMV')
          CALL ReadOptionR(buffer,1,opt_eir%dtimv) 
          opt_eir%dtimv = opt_eir%dtimv * 1.0E-06
        CASE('CON EIRENE TIME0')
          CALL ReadOptionR(buffer,1,opt_eir%time0) 
          opt_eir%time0 = opt_eir%time0 * 1.0E-06
        CASE('CON ITERATION BLOCK BEGIN')
          DO WHILE(osmGetLine(fp,buffer,ALL_LINES))            
c            WRITE(0,*) 'buffer 1 >'//TRIM(buffer)//'<'                        
            finished = .FALSE.
c...        Check if done:
            DO i = 1, LEN_TRIM(buffer)-25
c              WRITE(0,*) '>>'//buffer(i:i+24)//'<<'
              IF (buffer(i:i+24).EQ.'{CON ITERATION BLOCK END}') THEN
                finished = .TRUE.
                EXIT
              ENDIF
            ENDDO             
            IF (finished) EXIT
            niteration = niteration + 1
            IF (niteration.GT.MAX_NITERATION) 
     .        CALL ER('LoadControlOption','ITERATION_BUFFER overflow,'//
     .                ' input file likely malformed, sadly',*99)
            iteration_buffer(niteration) = buffer
          ENDDO
          IF (.NOT.finished) 
     .      CALL ER('LoadControlOption','Block end tag not found',*99)
        CASE DEFAULT
          CALL User_LoadOptions(fp,itag,buffer)
      ENDSELECT 

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE SplitBuffer(buffer,buffer_array)
      IMPLICIT none

      CHARACTER    , INTENT(IN)  :: buffer*(*)
      CHARACTER*256, INTENT(OUT) :: buffer_array(*)
c      CHARACTER, INTENT(OUT) :: buffer_array*256(*)  ! gfortran

      INTEGER i,j,k,n,m

      buffer_array(1) = ' '

      n = LEN_TRIM(buffer) + 1
      m = 0

      j = 0
      DO i = 1, n
        IF (buffer(i:i).EQ.' ' .AND.j.GT.0.OR.
     .      buffer(i:i).EQ.'"' .AND.j.LT.0.OR. 
     .      buffer(i:i).EQ.''''.AND.j.LT.0) THEN
          m = m + 1
c          WRITE(0,*) '>>>',j,i
          IF (j.LT.0) j = -j + 1
          buffer_array(m) = buffer(j:i-1)
c          WRITE(0,*) m,buffer(j:i-1)            
          j = 0
        ELSE
          IF (buffer(i:i).EQ.''''.AND.j.EQ.0) j = -i
          IF (buffer(i:i).EQ.'"' .AND.j.EQ.0) j = -i
          IF (buffer(i:i).NE.' ' .AND.j.EQ.0) j =  i
        ENDIF
      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadDivimpOption(fp,buffer,itag)
      USE mod_sol28_io
      USE mod_divimp
      IMPLICIT none

      INTEGER  , INTENT(IN)  :: fp,itag
      CHARACTER  :: buffer*(*)

      LOGICAL   osmGetLine

      INTEGER       i1,idum(5)
      LOGICAL       first_pass
      CHARACTER     cdum*1024  ! ,buffer_array*256(100) ! gfortran
      CHARACTER*256 buffer_array(100),buffer_list(100)
      REAL          version

      first_pass = .TRUE.

      DO i1 = 1, 100
        WRITE(buffer_array(i1),'(256X)')
      ENDDO

      SELECTCASE (buffer(3:itag-1))
c       ----------------------------------------------------------------
        CASE('DIV PARTICLE STATE')
          CALL ReadOptionI(buffer,1,opt_div%pstate)
c       ----------------------------------------------------------------
        CASE('DIV RIBBON GRID')
          READ(buffer(itag+2:itag+4),*) version
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
c            WRITE(0,*) 'buffer:'//TRIM(buffer)//'<'
            opt_div%rib_n = 1
            CALL SplitBuffer(buffer,buffer_array) 
            READ(buffer_array(1),*) opt_div%rib_type
            SELECTCASE (opt_div%rib_type)
              CASE(1)
                READ(buffer_array( 2),*) opt_div%rib_format
                opt_div%rib_file = TRIM(buffer_array(3))
              CASE(2)
                READ(buffer_array(2),*) opt_div%rib_rad_opt
                READ(buffer_array(3),*) opt_div%rib_rad_a
                READ(buffer_array(4),*) opt_div%rib_rad_d
                READ(buffer_array(5),*) opt_div%rib_rad_b
                READ(buffer_array(6),*) opt_div%rib_rad_c
              CASE(3)
                READ(buffer_array(2),*) opt_div%rib_pol_opt
                READ(buffer_array(3),*) opt_div%rib_pol_n
                READ(buffer_array(4),*) opt_div%rib_pol_a
                READ(buffer_array(5),*) opt_div%rib_pol_b
                READ(buffer_array(6),*) opt_div%rib_pol_c
                READ(buffer_array(7),*) opt_div%rib_pol_d
              CASE(4)
                READ(buffer_array(2),*) opt_div%rib_region
                READ(buffer_array(3),*) opt_div%rib_r1
                READ(buffer_array(4),*) opt_div%rib_r2
                READ(buffer_array(5),*) opt_div%rib_z1
                READ(buffer_array(6),*) opt_div%rib_z2
              CASE(5)  ! defaults
                READ(buffer_array(2),*) opt_div%rib_pol_n_def
                READ(buffer_array(3),*) opt_div%rib_pol_a_def
                READ(buffer_array(4),*) opt_div%rib_pol_b_def
                READ(buffer_array(5),*) opt_div%rib_pol_c_def
                READ(buffer_array(6),*) opt_div%rib_pol_d_def
              CASE DEFAULT
                CALL ER('LoadDivimpOption','Unknown ribbon grid '//
     .                  'option',*99)
            ENDSELECT
          ENDDO
c       ----------------------------------------------------------------
        CASE('DIV SPUTTER COMPILE')
          READ(buffer(itag+2:itag+4),*) version
          sputter_ndata = 0
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
c            WRITE(0,*) 'buffer:'//TRIM(buffer)//'<'
            sputter_ndata = sputter_ndata + 1
            buffer_list(sputter_ndata) = TRIM(buffer)
          ENDDO

          IF (ALLOCATED(sputter_data)) DEALLOCATE(sputter_data)
          ALLOCATE(sputter_data(sputter_ndata))

          DO i1 = 1, sputter_ndata
            buffer = TRIM(buffer_list(i1))
c           write(0,*) 'buffer: '//TRIM(buffer)
            CALL SplitBuffer(buffer,buffer_array) 
            READ(buffer_array(1),*) sputter_data(i1)%data_type
            SELECTCASE (sputter_data(i1)%data_type)
c             ----------------------------------------------------------
              CASE(1:3,5)
                sputter_data(i1)%case_name = TRIM(buffer_array(2))
                sputter_data(i1)%extension = TRIM(buffer_array(3))
                READ(buffer_array(4),*) sputter_data(i1)%fraction
                sputter_data(i1)%tag       = TRIM(buffer_array(5))
c             ----------------------------------------------------------
              CASE(4) ! sputtering specices is a constant fraction of the hydrogenic flux
                READ(buffer_array(2),*) sputter_data(i1)%atomic_number
                READ(buffer_array(3),*) sputter_data(i1)%atomic_mass
                READ(buffer_array(4),*) sputter_data(i1)%charge
                READ(buffer_array(5),*) sputter_data(i1)%fraction
                sputter_data(i1)%tag = TRIM(buffer_array(6))
                IF (sputter_data(i1)%fraction.EQ.-1.0) 
     .            sputter_data(i1)%fraction = 100.0
c             ----------------------------------------------------------
              CASE DEFAULT
c             ----------------------------------------------------------
                CALL ER('LoadDivimpOption','Unknown sputter data '//
     .                  'type',*99)
            ENDSELECT
          ENDDO          
c       ----------------------------------------------------------------
        CASE DEFAULT
          CALL User_LoadOptions(fp,itag,buffer)
      ENDSELECT 

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadEireneOption(fp,buffer,itag)
      USE mod_sol28_params
      USE mod_sol28_io
      USE mod_sol28_global
      USE mod_legacy
      USE mod_eirene06
      IMPLICIT none

      INTEGER  , INTENT(IN)  :: fp,itag
      CHARACTER  :: buffer*(*)

      LOGICAL   osmGetLine

      INTEGER   i1,idum(5)
      LOGICAL   first_pass
      CHARACTER     cdum*1024  ! ,buffer_array*256(100) ! gfortran
      CHARACTER*256 buffer_array(100)
      REAL      stratum_type,version,rdum(7),vec(3)

      first_pass = .TRUE.

      DO i1 = 1, 100
        WRITE(buffer_array(i1),'(256X)')
      ENDDO

      SELECTCASE (buffer(3:itag-1))
c       ----------------------------------------------------------------
        CASE('EIR IMPURITY SPUTTERING')
          CALL ReadOptionI(buffer,1,opt_eir%ilspt) 
c       ----------------------------------------------------------------
        CASE('EIR TET DUMP')
          CALL ReadOptionI(buffer,1,opt_eir%tet_iliin) 
c       ----------------------------------------------------------------
        CASE('EIR WHIPE')
          CALL ReadOptionI(buffer,1,opt_eir%whipe) 
c       ----------------------------------------------------------------
        CASE('EIR VOID GRID')
          opt_eir%nvoid = 0
          READ(buffer(itag+2:itag+4),*) opt_eir%void_version
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
c            WRITE(0,*) 'buffer:'//TRIM(buffer)//'<'
            opt_eir%nvoid = opt_eir%nvoid + 1
            IF     (opt_eir%void_version.EQ.1.0) THEN 
              READ(buffer,*) idum(1)
              IF (idum(1).EQ.-2) THEN
                READ(buffer,*) 
     .            opt_eir%void_zone(    opt_eir%nvoid),
     .            opt_eir%void_grid(1:2,opt_eir%nvoid)
              ELSE
                READ(buffer,*) 
     .            opt_eir%void_zone(    opt_eir%nvoid),
     .            opt_eir%void_grid(1:2,opt_eir%nvoid),
     .            opt_eir%void_wall(1:2,opt_eir%nvoid),
     .            opt_eir%void_add (1:2,opt_eir%nvoid),
     .            opt_eir%void_res (    opt_eir%nvoid),
     .            opt_eir%void_hole(1:2,opt_eir%nvoid),
     .            opt_eir%void_code(    opt_eir%nvoid),
     .            opt_eir%void_ne  (    opt_eir%nvoid),
     .            opt_eir%void_te  (    opt_eir%nvoid),
     .            opt_eir%void_ti  (    opt_eir%nvoid)
              ENDIF
            ELSEIF (opt_eir%void_version.EQ.2.0) THEN 
              CALL SplitBuffer(buffer,buffer_array) 
              i1 = opt_eir%nvoid
              READ(buffer_array(1),*) opt_eir%void_zone(i1)
              opt_eir%void2_grid(i1) = TRIM(buffer_array(2))
              opt_eir%void2_wall(i1) = TRIM(buffer_array(3))
              opt_eir%void2_add (i1) = TRIM(buffer_array(4))
              READ(buffer_array(5),*) opt_eir%void_res (i1)
              IF (LEN_TRIM(buffer_array(6)).GT.0) THEN
                READ(buffer_array(6),*) opt_eir%void_code(i1)
                READ(buffer_array(7),*) opt_eir%void_ne  (i1)
                READ(buffer_array(8),*) opt_eir%void_te  (i1)
                READ(buffer_array(9),*) opt_eir%void_ti  (i1)
              ENDIF
            ELSE
              CALL ER('LoadEireneOption','Unrecognised void '//
     .                'block version',*99)
            ENDIF            
          ENDDO
c       ----------------------------------------------------------------
        CASE('EIR ADDITIONAL SURFACES')
          opt_eir%nadd = 0
          READ(buffer(itag+2:itag+4),*) opt_eir%add_version
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
c            WRITE(0,*) 'buffer:'//TRIM(buffer)//'<'
            opt_eir%nadd = opt_eir%nadd + 1
            IF     (opt_eir%add_version.EQ.1.0) THEN 
              CALL SplitBuffer(buffer,buffer_array) 
c              write(0,*) '>'//TRIM(buffer_array(1))//'<'
c              write(0,*) '>'//TRIM(buffer_array(2))//'<'
c              write(0,*) '>'//TRIM(buffer_array(3))//'<'
c              write(0,*) '>'//TRIM(buffer_array(4))//'<'
c              write(0,*) '>'//TRIM(buffer_array(5))//'<'
              i1 = opt_eir%nadd
              READ(buffer_array(1),*) opt_eir%add_type (i1)
              READ(buffer_array(2),*) opt_eir%add_index(i1)
              SELECTCASE (opt_eir%add_type(i1))
                CASE(1)
                  opt_eir%add_file    (i1) = TRIM(buffer_array(3))
                  opt_eir%add_file_tag(i1) = TRIM(buffer_array(4))
                  opt_eir%add_tag     (i1) = TRIM(buffer_array(5))
                CASE(2)
                  READ(buffer_array( 3),*) opt_eir%add_holex(i1)
                  READ(buffer_array( 4),*) opt_eir%add_holey(i1)
                  opt_eir%add_tag     (i1) = TRIM(buffer_array(5))
                CASE DEFAULT
                  CALL ER('LoadEireneOption','Unknown additional '//
     .                    'surface type',*99)
              ENDSELECT
            ELSE
              CALL ER('LoadEireneOption','Unrecognised add '//
     .                'block version',*99)
            ENDIF            
          ENDDO
c       ----------------------------------------------------------------
        CASE('E NEUTRAL SOURCES')
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
            IF (first_pass) THEN
              opt_eir%nstrata = 0
              first_pass = .FALSE.
            ENDIF
            opt_eir%nstrata = opt_eir%nstrata + 1
c            WRITE(0,*) 'BUFFER:',TRIM(buffer)
            READ(buffer,*) 
     .        opt_eir%type         (opt_eir%nstrata),
     .        opt_eir%npts         (opt_eir%nstrata),
     .        opt_eir%flux         (opt_eir%nstrata),
     .        opt_eir%flux_fraction(opt_eir%nstrata),
     .        opt_eir%species      (opt_eir%nstrata),
     .        opt_eir%species_index(opt_eir%nstrata),
     .        opt_eir%sorene       (opt_eir%nstrata)
            IF     (opt_eir%type(opt_eir%nstrata).EQ.1.0) THEN  ! Target surface flux
              READ(buffer,*) rdum(1:7),
     .          opt_eir%target(opt_eir%nstrata),
     .          opt_eir%txtsou(opt_eir%nstrata)
              opt_eir%range_tube(1,opt_eir%nstrata) = 1
              opt_eir%range_tube(2,opt_eir%nstrata) = 99999
            ELSEIF (opt_eir%type(opt_eir%nstrata).EQ.1.1) THEN  ! Target surface flux
              READ(buffer,*) rdum(1:7),
     .          opt_eir%target        (opt_eir%nstrata),
     .          opt_eir%range_tube(1:2,opt_eir%nstrata),
     .          opt_eir%txtsou        (opt_eir%nstrata)
            ELSEIF (opt_eir%type(opt_eir%nstrata).EQ.2.0) THEN  ! Volume recombination
              READ(buffer,*) rdum(1:7),
     .          opt_eir%txtsou(opt_eir%nstrata)
            ELSEIF (opt_eir%type(opt_eir%nstrata).EQ.3.0.OR.
     .              opt_eir%type(opt_eir%nstrata).EQ.3.1.OR.    ! Point source injection (gas puff, beams)
     .              opt_eir%type(opt_eir%nstrata).EQ.4.0) THEN  ! Surface
              READ(buffer,*) rdum(1:7),
     .          opt_eir%sorcos   (opt_eir%nstrata),
     .          opt_eir%sormax   (opt_eir%nstrata),
     .          opt_eir%sorad(1:6,opt_eir%nstrata),
     .          opt_eir%txtsou   (opt_eir%nstrata)

              IF (opt_eir%sorad(4,opt_eir%nstrata).EQ.0.0.AND.  ! Make sure that the launch vector is not purely
     .            opt_eir%sorad(5,opt_eir%nstrata).NE.0.0.AND.  ! vertical, since the 202 launch distribution
     .            opt_eir%sorad(6,opt_eir%nstrata).EQ.0.0)      ! is hardwired in SetupEireneStrata
     .          opt_eir%sorad(4,opt_eir%nstrata) = 1.0E-06 

            ELSE
              CALL ER('LoadEireneOption','Unknown stratum type',*99)
            ENDIF
          ENDDO            
c          WRITE(0,*) 'OPT_EIR%NSTRATA:',opt_eir%strata,rdum(1:6)
c          WRITE(0,*) 'OPT_EIR%NSTRATA:',opt_eir%nstrata,
c     .                              opt_eir%sorad(opt_eir%nstrata)
c          STOP
c       ----------------------------------------------------------------
        CASE('EIR PARTICLE SPECTRA')
          opt_eir%nadspc = 0
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
            opt_eir%nadspc = opt_eir%nadspc + 1
            READ(buffer,*) 
     .        opt_eir%isrfcll   (opt_eir%nadspc),  ! Type of spectrum, 0=surface flux, 1=cell average                    
     .        opt_eir%ispsrf    (opt_eir%nadspc),  ! Surface/cell index, <0=non-default standard, >0=additional surfaces          
     .        opt_eir%ispsrf_ref(opt_eir%nadspc),  ! Which code does the surface index refer to?                             
     .        opt_eir%iptyp     (opt_eir%nadspc),  ! Species type eg 1=atoms, 2=molecules, 3=test ions, 4=?                  
     .        opt_eir%ipsp      (opt_eir%nadspc),  ! Species sub-index eg, 1=first atom species, 2=second atom species, etc.   
     .        opt_eir%isptyp    (opt_eir%nadspc),  ! Spectrum type wrt units, 1=1/eV/s, 2=1/s                                  
     .        opt_eir%nsps      (opt_eir%nadspc),  ! Number of bins                                                          
     .        opt_eir%spcmn     (opt_eir%nadspc),  ! Lower bound of energy range for spectrum                                
     .        opt_eir%spcmx     (opt_eir%nadspc),  ! Upper bound                                                             
     .        opt_eir%idirec    (opt_eir%nadspc)   ! If >0 then a projection on a direction is used in the statistics (??)   
            SELECTCASE (opt_eir%idirec(opt_eir%nadspc)) ! 1=vector for projecting onto, 2=collect cells along a LOS, 3=same, but project onto vector as well
              CASE (0)
              CASE (1)
                READ(buffer,*) idum(1:2),cdum   ,idum(1:4),
     .                         rdum(1:2),idum(1),
     .            opt_eir%spcvx(opt_eir%nadspc),     ! Don't know really, but it was in the example that AK sent, originally from VK
     .            opt_eir%spcvy(opt_eir%nadspc),
     .            opt_eir%spcvz(opt_eir%nadspc)
              CASE (2:3)                             
                READ(buffer,*) idum(1:2),cdum   ,idum(1:4),
     .                         rdum(1:2),idum(1),
     .            opt_eir%spc_p1(opt_eir%nadspc,1:3),  ! Vector along which volume (cell) energy distributions are requested
     .            opt_eir%spc_p2(opt_eir%nadspc,1:3)

                IF (opt_eir%idirec(opt_eir%nadspc).EQ.3) THEN
                  vec(1:3) = opt_eir%spc_p2(opt_eir%nadspc,1:3) -
     .                       opt_eir%spc_p1(opt_eir%nadspc,1:3)
                  vec = vec / SQRT(vec(1)**2+vec(2)**2+vec(3)**2)
                  opt_eir%spcvx(opt_eir%nadspc) = vec(1) 
                  opt_eir%spcvy(opt_eir%nadspc) = vec(2)
                  opt_eir%spcvz(opt_eir%nadspc) = vec(3)
                ENDIF
              CASE DEFAULT
                CALL ER('LoadEireneOption','Unknown spectrum '//
     .                  'direction',*99)
            ENDSELECT
          ENDDO
c       ----------------------------------------------------------------
        CASE('EIR TETRAHEDRON GRID')
          opt_eir%tet_n = 0
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
            opt_eir%tet_n = opt_eir%tet_n + 1
            i1 = opt_eir%tet_n
            CALL SplitBuffer(buffer,buffer_array) 
            SELECTCASE (TRIM(buffer_array(1)))
              CASE ('1.0')
                READ(buffer,*) 
     .            opt_eir%tet_type(i1),  !                
     .            opt_eir%tet_x1  (i1),  ! 
     .            opt_eir%tet_y1  (i1),  ! 
     .            opt_eir%tet_x2  (i1),  ! 
     .            opt_eir%tet_y2  (i1)   ! 
              CASE ('2.0')
                READ(buffer_array(1),*) opt_eir%tet_type    (i1)  !         
                READ(buffer_array(2),*) opt_eir%tet_index   (i1)  !         
                READ(buffer_array(3),*) opt_eir%tet_mode    (i1)  !         
                READ(buffer_array(4),*) opt_eir%tet_param1  (i1)  !         
                READ(buffer_array(5),*) opt_eir%tet_param2  (i1)  !         
                opt_eir%tet_del_hole(i1) = TRIM(buffer_array(6))  !         
                opt_eir%tet_del_zone(i1) = TRIM(buffer_array(7))  !         
              CASE ('3.0')
                READ(buffer_array(1),*) opt_eir%tet_type    (i1)  !         
                READ(buffer_array(2),*) opt_eir%tet_index   (i1)  !         
                opt_eir%tet_sec_list(i1) = TRIM(buffer_array(3))  !         
              CASE ('4.0')
                READ(buffer_array(1),*) opt_eir%tet_type     (i1)  !         
                opt_eir%tet_composite(i1) = TRIM(buffer_array(2))  !         
                READ(buffer_array(3),*) opt_eir%tet_offset   (i1)  !         
              CASE ('5.0')  ! Grid refinement
                READ(buffer_array(1),*) opt_eir%tet_type(i1)  !         
                READ(buffer_array(2),*) opt_eir%tet_mode(i1)  !         
                SELECTCASE (opt_eir%tet_mode(i1))
                  CASE (4)
                    READ(buffer_array(3 ),*) opt_eir%tet_param1(i1)  !         
                    READ(buffer_array(4 ),*) opt_eir%tet_param2(i1)  !         
                    READ(buffer_array(5 ),*) opt_eir%tet_param3(i1)  !         
                    READ(buffer_array(6 ),*) opt_eir%tet_x1    (i1)  !         
                    READ(buffer_array(7 ),*) opt_eir%tet_y1    (i1)  !         
                    READ(buffer_array(8 ),*) opt_eir%tet_z1    (i1)  !         
                    READ(buffer_array(9 ),*) opt_eir%tet_x2    (i1)  !         
                    READ(buffer_array(10),*) opt_eir%tet_y2    (i1)  !         
                    READ(buffer_array(11),*) opt_eir%tet_z2    (i1)  !         
                  CASE DEFAULT
                    CALL ER('LoadEireneOption','Madness, no idea '//
     .                      'what''s going on',*99)
                ENDSELECT
              CASE DEFAULT
                CALL ER('LoadEireneOption','Unknown tetrahedron grid '//
     .                  'TYPE found',*99)
            ENDSELECT
          ENDDO
c       ----------------------------------------------------------------
        CASE('EIR SURFACE PROPERTIES')
          opt_eir%sur_n  = 0
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
c            WRITE(0,*) 'buffer:'//TRIM(buffer)
            opt_eir%sur_n = opt_eir%sur_n + 1
            CALL SplitBuffer(buffer,buffer_array) 
            READ(buffer_array( 1),*) opt_eir%sur_type  (opt_eir%sur_n)
            opt_eir%sur_index (opt_eir%sur_n) = TRIM(buffer_array(2))
            opt_eir%sur_sector(opt_eir%sur_n) = TRIM(buffer_array(3))
            READ(buffer_array( 4),*) opt_eir%sur_iliin (opt_eir%sur_n)
            READ(buffer_array( 5),*) opt_eir%sur_ilside(opt_eir%sur_n)
            READ(buffer_array( 6),*) opt_eir%sur_ilswch(opt_eir%sur_n)
            READ(buffer_array( 7),*) opt_eir%sur_tr1   (opt_eir%sur_n)
            READ(buffer_array( 8),*) opt_eir%sur_tr2   (opt_eir%sur_n)
            READ(buffer_array( 9),*) opt_eir%sur_recyct(opt_eir%sur_n)
            READ(buffer_array(10),*) opt_eir%sur_ilspt (opt_eir%sur_n)
            READ(buffer_array(11),*) opt_eir%sur_temp  (opt_eir%sur_n)
            opt_eir%sur_mat(opt_eir%sur_n) = TRIM(buffer_array(12))
            READ(buffer_array(13),*) opt_eir%sur_hard  (opt_eir%sur_n)
            READ(buffer_array(14),*) opt_eir%sur_remap (opt_eir%sur_n)
            opt_eir%sur_tag(opt_eir%sur_n) = TRIM(buffer_array(15))
          ENDDO
c       ----------------------------------------------------------------
       CASE DEFAULT
          CALL User_LoadOptions(fp,itag,buffer)
      ENDSELECT 

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadFilamentOption(fp,buffer,itag)
      USE mod_sol28_params
      USE mod_sol28_io
      USE mod_options
      USE mod_legacy
      USE mod_eirene06
      IMPLICIT none

      INTEGER  , INTENT(IN)  :: fp,itag
      CHARACTER  :: buffer*(*)

      LOGICAL   osmGetLine

      INTEGER   i1
      REAL      stratum_type,version,rdum(7),length(2)

      SELECTCASE (buffer(3:itag-1))
        CASE('FIL OPT')
          CALL ReadOptionI(buffer,1,opt_fil%opt)
        CASE('FIL TARGET FLUX')
          CALL ReadOptionI(buffer,1,opt_fil%target_flux)
        CASE('FIL LENGTH')
c          WRITE(0,*) 'BUFFER:',TRIM(buffer)
          CALL ReadOptionR(buffer,2,length)
          opt_fil%length1 = length(1)
          opt_fil%length2 = length(2)
        CASE('FIL CLIP')
          CALL ReadOptionI(buffer,1,opt_fil%clip)
        CASE('FIL TIME OF INTEREST','FIL START TIME')
          CALL ReadOptionR(buffer,1,opt_fil%start_time)
          opt_fil%start_time = opt_fil%start_time * 1.0E-06
        CASE('FIL TIME STEP')
          CALL ReadOptionR(buffer,1,opt_fil%time_step)
          opt_fil%time_step = opt_fil%time_step * 1.0E-06
        CASE DEFAULT
          CALL User_LoadOptions(fp,itag,buffer)
      ENDSELECT 

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadFileOption(fp,buffer,itag)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      INTEGER  , INTENT(IN) :: fp,itag
      CHARACTER, INTENT(IN) :: buffer*(*)
 
      INTEGER   i
      CHARACTER cdum1*512

      SELECTCASE (buffer(3:itag-1))
        CASE('F LOG')
          READ(buffer,*) cdum1,logop
          IF (logop.LT.0) THEN
            CLOSE(logfp)
            logfp = 0
            logop = -logop
          ENDIF
          opt%log   = logop
          opt%logfp = logfp
        CASE('F OSM_LOAD')
          opt%osm_load = 1
          READ(buffer,*) cdum1,opt%f_osm_load
        CASE('F OSM_DIR')
          READ(buffer,*) cdum1,opt%f_osm_dir
        CASE('F EIRENE_DIR')
          READ(buffer,*) cdum1,opt_eir%f_eirene_dir
        CASE('F EIRENE_13')
          opt_eir%f_eirene_load = 1
          READ(buffer,*) cdum1,opt_eir%f_eirene_13
        CASE('F EIRENE_15')
          opt_eir%f_eirene_load = 1
          READ(buffer,*) cdum1,opt_eir%f_eirene_15
c          WRITE(0,*) 'opt%f_eirene_15:',TRIM(opt_eir%f_eirene_15)
        CASE('GRID FORMAT')
          CALL ReadOptionI(buffer,1,opt%f_grid_format)
        CASE('GRID LOAD METHOD')
          CALL ReadOptionI(buffer,1,opt%f_grid_load_method)
        CASE('GRID FILE')
          READ(buffer,*) cdum1,opt%f_grid_file
        CASE('GRID STRIP CELLS')
          CALL ReadOptionI(buffer,1,opt%f_grid_strip)
        CASE('GRID DEL TUBES')
          READ(buffer,*,END=10) cdum1,opt%grd_tdel(1:S28_MAXNTDEL)
 10       DO i = 1, S28_MAXNTDEL
            IF (opt%grd_tdel(i).NE.0) opt%grd_ntdel = opt%grd_ntdel + 1
          ENDDO
c          WRITE(0,*) 'GRID :',opt%grd_tdel(1:opt%grd_ntdel)
        CASE('FILE MATERIALS DATA')
          READ(buffer,*) cdum1,opt%mat_opt,opt%mat_file
      CASE DEFAULT
          CALL User_LoadOptions(fp,itag,buffer)
      ENDSELECT 

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadSourceOption(fp,buffer,itag)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      INTEGER  , INTENT(IN) :: fp,itag
      CHARACTER, INTENT(IN) :: buffer*(*)

      SELECTCASE (buffer(3:itag-1))
        CASE('S P_ION')
          CALL ReadOptionI(buffer,2,opt%p_ion)
        CASE('S P_ION_EXP')
          CALL ReadOptionR(buffer,2,opt%p_ion_exp)
        CASE('S P_ION_FRAC')
          CALL ReadOptionR(buffer,2,opt%p_ion_frac)
        CASE('S P_REC')
          CALL ReadOptionI(buffer,2,opt%p_rec)
        CASE('S P_ANO')
          CALL ReadOptionI(buffer,2,opt%p_ano)
        CASE('S P_ANO_DIST')
          CALL ReadOptionI(buffer,2,opt%p_ano_dist)
        CASE('S P_ANO_EXP')
          CALL ReadOptionR(buffer,2,opt%p_ano_exp)
        CASE('S M_PIN')
          CALL ReadOptionI(buffer,2,opt%m_mom)
        CASE('S M_FIT')
          CALL ReadOptionI(buffer,2,opt%m_fit)
        CASE('S M_ANO')
          CALL ReadOptionI(buffer,2,opt%m_ano)
        CASE('S M_ANO_DIST')
          CALL ReadOptionI(buffer,2,opt%m_ano_dist)
        CASE('S M_ANO_EXP')
          CALL ReadOptionR(buffer,2,opt%m_ano_exp)
        CASE('S TE_ION')
          CALL ReadOptionI(buffer,2,opt%te_ion)
        CASE('S TE_REC')
          CALL ReadOptionI(buffer,2,opt%te_rec)
        CASE('S TE_ANO')
          CALL ReadOptionI(buffer,2,opt%te_ano)
        CASE('S TE_ANO_PSOL')
          CALL ReadOptionI(buffer,2,opt%te_ano_psol)
        CASE('S TI_ION')
          CALL ReadOptionI(buffer,2,opt%ti_ion)
        CASE('S TI_REC')
          CALL ReadOptionI(buffer,2,opt%ti_rec)
        CASE('S TI_ANO')
          CALL ReadOptionI(buffer,2,opt%ti_ano)
      CASE DEFAULT
          CALL User_LoadOptions(fp,itag,buffer)
      ENDSELECT 

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadSolverOption(fp,buffer,itag)
      USE mod_sol28_params
      USE mod_sol28_io
      USE mod_sol28_global
      IMPLICIT none

      INTEGER  , INTENT(IN) :: fp,itag
      CHARACTER, INTENT(IN) :: buffer*(*)

      LOGICAL osmGetLine

      CHARACTER*256 buffer_array(100)  ! gfortran

      SELECTCASE (buffer(3:itag-1))
c       ----------------------------------------------------------------
        CASE('SOL APPLICATION')
          opt%sol_n = 0
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
            opt%sol_n = opt%sol_n + 1
            READ(buffer,*) opt%sol_tube  (1:2,opt%sol_n),
     .                     opt%sol_option(    opt%sol_n)
          ENDDO
c       ----------------------------------------------------------------
        CASE('SOL RADIAL VELOCITY')
          CALL SplitBuffer(buffer,buffer_array) 
          READ(buffer_array(2),*) opt%radvel
          SELECTCASE (opt%radvel)
            CASE (0) 
            CASE (1) 
              READ(buffer_array(3),*) opt%radvel_param(1)
            CASE DEFAULT
              CALL ER('LoadSolverOption','Invalid radial velocity, '//
     .                'option',*99)
          ENDSELECT
c       ----------------------------------------------------------------
        CASE DEFAULT
          CALL User_LoadOptions(fp,itag,buffer)
c       ----------------------------------------------------------------
      ENDSELECT 

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadSOLPSOption(fp,buffer,itag)
      USE mod_sol28_io
      USE mod_solps_params
      USE mod_solps
      IMPLICIT none

      INTEGER  , INTENT(IN) :: fp,itag
      CHARACTER, INTENT(IN) :: buffer*(*)

      CHARACTER cdum1*1024

      LOGICAL osmGetLine

      SELECTCASE (buffer(3:itag-1))
        CASE('SOLPS LOAD SOLUTION')
          nsolps_data = 0
          IF (ALLOCATED(solps_data)) THEN
            CALL WN('LoadSOLPSOptions','SOLPS data array '//
     .              'already allocated, deallocating')
            DEALLOCATE(solps_data)
          ENDIF
          READ(buffer,*) cdum1,solps_opt
          ALLOCATE(solps_data(MAX_SOLPS_DATA))
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
            nsolps_data = nsolps_data + 1
            IF (nsolps_data.GT.MAX_SOLPS_DATA) 
     .        CALL ER('LoadSOLPSOption','Too many data lines, '//
     .                'increase MAX_SOLPS_DATA in mod_solps.f',*99)
            READ(buffer,*) solps_data(nsolps_data)%fname,
     .                     solps_data(nsolps_data)%format,
     .                     solps_data(nsolps_data)%column,
     .                     solps_data(nsolps_data)%type,
     .                     solps_data(nsolps_data)%tag,
     .                     solps_data(nsolps_data)%z,
     .                     solps_data(nsolps_data)%a,
     .                     solps_data(nsolps_data)%charge
          ENDDO
      CASE DEFAULT
        CALL User_LoadOptions(fp,itag,buffer)
      ENDSELECT 

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadTransportOption(fp,buffer,itag)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      INTEGER  , INTENT(IN) :: fp,itag
      CHARACTER, INTENT(IN) :: buffer*(*)

      SELECTCASE (buffer(3:itag-1))
        CASE('T TE_KAPPA')
          CALL ReadOptionR(buffer,2,opt%te_kappa)
        CASE('T TE_CONV')
          CALL ReadOptionI(buffer,2,opt%te_conv)
        CASE('T TE_FLUXLIMIT')
          CALL ReadOptionI(buffer,2,opt%te_fluxlimit)
        CASE('T TI')
          CALL ReadOptionI(buffer,2,opt%ti)
        CASE('T TI_KAPPA')
          CALL ReadOptionR(buffer,2,opt%ti_kappa)
        CASE('T TI_CONV')
          CALL ReadOptionI(buffer,2,opt%ti_conv)
        CASE('T TI_RATIO')
          CALL ReadOptionR(buffer,2,opt%ti_ratio)
        CASE('T TI_EQUIL')
          CALL ReadOptionI(buffer,2,opt%ti_equil)
        CASE('T BC')
          CALL ReadOptionI(buffer,2,opt%bc)
        CASE('T SUPER')
          CALL ReadOptionI(buffer,2,opt%super)
      CASE DEFAULT
          CALL User_LoadOptions(fp,itag,buffer)
      ENDSELECT 

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadWallOption(fp,buffer,itag)
      USE mod_sol28_io
      USE mod_sol28_wall
      IMPLICIT none

      INTEGER  , INTENT(IN) :: fp,itag
      CHARACTER, INTENT(IN) :: buffer*(*)

      CHARACTER cdum1*1024

      LOGICAL osmGetLine

      SELECTCASE (buffer(3:itag-1))
        CASE('WALL 2D SEGMENTS')
          nopt_wall = 0
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
            nopt_wall = nopt_wall + 1
            READ(buffer,*) opt_wall(nopt_wall)%type,
     .                     opt_wall(nopt_wall)%class,
     .                     opt_wall(nopt_wall)%index(WAL_GROUP),
     .                     opt_wall(nopt_wall)%material_tag,
     .                     opt_wall(nopt_wall)%temperature,
     .                     opt_wall(nopt_wall)%index(WAL_RANGE1),
     .                     opt_wall(nopt_wall)%index(WAL_RANGE2),
     .                     opt_wall(nopt_wall)%file_format,
     .                     opt_wall(nopt_wall)%file_name
          ENDDO
      CASE DEFAULT
        CALL User_LoadOptions(fp,itag,buffer)
      ENDSELECT 

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadMiscOption(fp,buffer,itag,status)
      USE mod_sol28_params
      USE mod_sol28_io
      USE mod_sol28_global
      USE mod_legacy
      IMPLICIT none

      INTEGER  , INTENT(IN)  :: fp,itag
      CHARACTER  :: buffer*(*)
c      CHARACTER, INTENT(IN)  :: buffer*(*)
      LOGICAL  , INTENT(OUT) :: status

      LOGICAL osmGetLine

      INTEGER   i1,sub_option
      LOGICAL   node_fit,node_data,ldum1
      REAL      node_type,rdum1,rdum2
      CHARACTER cdum1
      TYPE(type_node) :: node_tmp

      status = .TRUE.

      SELECTCASE (buffer(3:itag-1))
        CASE('MATERIAL DATA')

        CASE('088')
          tarninter(HI) = 0
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
            tarninter(HI) = tarninter(HI) + 1
            CALL PadBufferR(15,buffer,0.0)
            READ(buffer,*) tarinter(tarninter(HI),1:4,HI)
          ENDDO
        CASE('089')
          tarninter(LO) = 0
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
            tarninter(LO) = tarninter(LO) + 1
            CALL PadBufferR(15,buffer,0.0)
            READ(buffer,*) tarinter(tarninter(LO),1:4,LO)
          ENDDO
        CASE('S74')
          node_data = .FALSE.
          READ(buffer(itag+2:itag+4),*) opt%s28mode  
          IF    (opt%s28mode.EQ.4.0) THEN
            osmns28 = 0            
            DO WHILE(osmGetLine(fp,buffer,NO_TAG))
              osmns28 = osmns28 + 1
              READ(buffer,*) osms28(osmns28,1:12)
            ENDDO
          ELSEIF (opt%s28mode.EQ.4.1) THEN
            node_fit = .FALSE.
            osmnnode = 0            
            DO WHILE(osmGetLine(fp,buffer,NO_TAG))
              CALL PadBufferR(15,buffer,0.0)
              READ(buffer,*) node_type
              osmnnode = osmnnode + 1
              IF (node_type.GE.1.0.AND.node_type.LE.4.0) THEN
c...            Load main node interpolation data:
                osmnode(osmnnode)%type = 0.0  ! Add spacer to separate data cards
                READ(buffer,*) 
     .            node_tmp%type,
     .            node_tmp%tube_range(1:2),
     .            node_tmp%rad_mode,
     .            node_tmp%rad_coord,
     .            node_tmp%rad_exp,
     .            node_tmp%par_mode,
     .            node_tmp%par_exp,
     .            node_tmp%par_set
                IF (node_tmp%rad_mode.EQ.4) THEN
                  ldum1 = osmGetLine(fp,buffer,NO_TAG)
                  CALL PadBufferR(15,buffer,1.0)
                  READ(buffer,*) node_tmp%file_name,
     .                           node_tmp%file_format,
     .                           node_tmp%file_shift,
     .                           node_tmp%file_scale_ne,
     .                           node_tmp%file_scale_M,
     .                           node_tmp%file_scale_pe,
     .                           node_tmp%file_scale_Te,
     .                           node_tmp%file_scale_Ti,
     .                           node_tmp%file_scale_V
                ENDIF
                node_data = .TRUE.
                node_fit  = .FALSE.
c                IF (node_tmp%rad_mode.EQ.5) THEN
c                  node_fit  = .TRUE.
c                  node_data = .FALSE.
c                ENDIF
              ELSEIF (node_type.EQ.0.0.AND.node_data) THEN
                osmnode(osmnnode) = node_tmp
                osmnode(osmnnode)%fit_type = 0.0
                READ(buffer,*) rdum1,
     .            osmnode(osmnnode)%rad_x,       ! Change to just x and y?
     .            osmnode(osmnnode)%rad_y,
     .            osmnode(osmnnode)%ne,
     .            osmnode(osmnnode)%v,
     .            osmnode(osmnnode)%pe,
     .            osmnode(osmnnode)%te,
     .            osmnode(osmnnode)%ti(1),       ! Assumption
     .            osmnode(osmnnode)%epot
              ELSEIF (node_type.EQ.-1.0.AND.node_data) THEN
                osmnode(osmnnode)%type = 0.0
                IF (node_tmp%rad_mode.EQ.6) THEN
c...              Load data for radial fits to pedestal prescription:
                  READ(buffer,*) rdum1,rdum2
                  SELECTCASE (NINT(rdum2))                 
                    CASE (-1) 
                    CASE ( 1)  ! Linear in core, TANH in pedestal, exponential in SOL
                      READ(buffer,*) rdum1,
     .                  osmnode(osmnnode)%fit_type,
     .                  osmnode(osmnnode)%fit_quantity

                      sub_option = NINT(10 * 
     .                  MOD(          osmnode(osmnnode)%fit_quantity,
     .                      REAL(NINT(osmnode(osmnnode)%fit_quantity))))
c                      write(0,*) 'sub_option',sub_option
                      SELECTCASE (sub_option)
                        CASE (1)  ! Adding adjustmend of Ti:Te ratio cross-over location
                          osmnode(osmnnode)%fit_width = 2
                          READ(buffer,*) rdum1,rdum1,rdum1,
     .                      osmnode(osmnnode)%fit_p(1:10)
                        CASE DEFAULT
                          osmnode(osmnnode)%fit_width = 1
                          READ(buffer,*) rdum1,rdum1,rdum1,
     .                      osmnode(osmnnode)%fit_p(1:9)
                          osmnode(osmnnode)%fit_p(10) = 3.0  ! default, for backward compatibility
                      ENDSELECT
c                      READ(buffer,*) rdum1,
c     .                  osmnode(osmnnode)%fit_type,
c     .                  osmnode(osmnnode)%fit_quantity,
c     .                  osmnode(osmnnode)%fit_p(1:9)

                    CASE ( 2)  ! Linear in the core, exponential in SOL
                      READ(buffer,*) rdum1,
     .                  osmnode(osmnnode)%fit_type,
     .                  osmnode(osmnnode)%fit_quantity,
     .                  osmnode(osmnnode)%fit_p(1:7)
                    CASE DEFAULT
                      CALL ER('LoadMiscOption','Bad TYPE',*99)
                  ENDSELECT
                ELSE
c...              Load data for radial fits to functions:
                  READ(buffer,*) rdum1,
     .              osmnode(osmnnode)%fit_type,
     .              osmnode(osmnnode)%fit_psin(1:2),
     .              osmnode(osmnnode)%fit_shift,
     .              osmnode(osmnnode)%fit_quantity,
     .              osmnode(osmnnode)%fit_p(1:6)
                ENDIF
              ELSEIF (node_type.EQ.0.0) THEN
c...            Spacer, ignore:
                osmnode(osmnnode)%type = 0.0
              ELSE
                CALL ER('LoadMiscOption','Unrecognized node type',*99)
              ENDIF
            ENDDO            
          ELSE
            CALL ER('LoadMiscOption','Unrecognized S28MODE value',*99)
          ENDIF


          WRITE(88,*) 'DEBUG 11:',osmnode(11)%type,osmnode(11)%fit_type


          DO i1 = 1, osmnnode
            WRITE(logfp,*) 'nodes:',osmnode(i1)%type,
     .                     osmnode(i1)%par_mode,
     .                     osmnode(i1)%v,
     .                     osmnode(i1)%file_scale_M
          ENDDO

        CASE('030')
        CASE('E16')
        CASE('999','EXIT')
          status = .FALSE.
       CASE DEFAULT
          CALL User_LoadOptions(fp,itag,buffer)
      ENDSELECT 

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE osm_InitializeOptions
      USE mod_sol28_global
      USE mod_options
      USE mod_eirene06
      USE mod_eirene_history
      USE mod_legacy
      USE mod_solps
      USE mod_divimp
      IMPLICIT none

      CALL InitializeLegacyVariables

c...  Input variables:
      iiteration = 0
      niteration = 0

c...  Output options:
      logop     = 1
      logfp     = 88
      opt%log   = logop
      opt%logfp = logfp

      opt%osm_load = 0

c...  Node interpolation data:
      opt%s28mode = 4.0
      osmnnode = 0

c...  

c...  OSM options:
      nion = 1

      opt%osm_load = 0
      WRITE(opt%f_osm_dir ,'(512X)')
      WRITE(opt%f_osm_load,'(512X)')
      opt_eir%f_eirene_load = 0
      WRITE(opt_eir%f_eirene_dir,'(512X)')
      WRITE(opt_eir%f_eirene_13 ,'(512X)')
      WRITE(opt_eir%f_eirene_15 ,'(512X)')

      opt%f_grid_format = 0
      opt%f_grid_load_method = 2 
      WRITE(opt%f_grid_file,'(512X)')
      opt%f_grid_strip = 0
      opt%grd_ntdel = 0
      opt%grd_tdel  = 0


      opt%mat_opt  = 1
      opt%mat_file = 'materials.dat'

      opt%sol_n = 1
      opt%sol_tube(1,1) = 1
      opt%sol_tube(2,1) = 1E+8
      opt%sol_option(1) = 28

      opt%pin_data   = .FALSE.
      opt%p_ion      = 2   ! changed from 3 on 23/07/2010
      opt%p_ion_exp  = 0.1
      opt%p_ion_frac = 100.0
      opt%p_rec      = 0 
      opt%p_ano      = 2
      opt%p_ano_dist = 3   ! new on 23/04/2012
      opt%p_ano_exp  = 0.0 ! new on 23/04/2012
      opt%m_mom      = 0
      opt%m_fit      = 2
      opt%m_ano      = 2
      opt%m_ano_dist = 3   ! changed from 1 on 23/07/2010
      opt%m_ano_exp  = 2.0 ! changed from 0.0 on 23/07/2010

      opt%te_rec = 0
      opt%te_ion = 0
      opt%te_ano = 0
      opt%te_ano_psol = 0
      opt%te_fluxlimit = 0
      opt%te_kappa = 2000.0
      opt%te_conv = 1

      opt%ti = 0
      opt%ti_rec = 0
      opt%ti_ion = 0
      opt%ti_ano = 0
      opt%ti_ano_psol = 0
      opt%ti_kappa = 60.0
      opt%ti_conv = 0
      opt%ti_ratio = 1.0
      opt%ti_equil = 0

      opt%bc    = 1
      opt%super = 0

      opt%radvel = 0

      ref_nion   = 0
      ref_nfluid = 1

c...  Filament options:
      opt_fil%opt = 0
      opt_fil%target_flux = 0
      opt_fil%clip = 0
      opt_fil%start_time = 0.0
      opt_fil%time_step = 10.0  ! (microseconds)
      opt_fil%scale(1) = 5.0
      opt_fil%scale(2) = 2.0
      opt_fil%scale(3) = 1.0
      opt_fil%length1 = -99.0
      opt_fil%length2 = -99.0

c...  Eirene options:

c      opt_eir%nvoid = 0
      opt_eir%nvoid = 1
      opt_eir%void_version   =  1.0
      opt_eir%void_zone(  1) =   -1
      opt_eir%void_grid(1,1) =    2
      opt_eir%void_grid(2,1) =  999
      opt_eir%void_wall(:,1) =   -1
      opt_eir%void_add (:,1) =   -1
      opt_eir%void_hole(:,1) = -1.0
      opt_eir%void_res (  1) =  0.1
      opt_eir%void_code(  1) =   -1
      opt_eir%void_ne  (  1) =  0.0
      opt_eir%void_te  (  1) =  0.0
      opt_eir%void_ti  (  1) =  0.0
      opt_eir%void_tag (  1) =  'default'

      opt_eir%nadd = 0

      opt_eir%nstrata = 3
      opt_eir%type         (1) = 1.0
      opt_eir%npts         (1) = -90000
      opt_eir%flux         (1) = 1.0
      opt_eir%flux_fraction(1) = 1.0
      opt_eir%species      (1) = 4
      opt_eir%species_index(1) = 1
      opt_eir%sorene       (1) = 0.0
      opt_eir%target       (1) = 1
      opt_eir%txtsou       (1) = 'default inner target'
      opt_eir%range_tube (1,1) = 1
      opt_eir%range_tube (2,1) = 99999

      opt_eir%type         (2) = 1.0
      opt_eir%npts         (2) = -90000
      opt_eir%flux         (2) = 1.0
      opt_eir%flux_fraction(2) = 1.0
      opt_eir%species      (2) = 4
      opt_eir%species_index(2) = 1
      opt_eir%sorene       (2) = 0.0
      opt_eir%target       (2) = 2
      opt_eir%txtsou       (2) = 'default outer target'
      opt_eir%range_tube (1,2) = 1
      opt_eir%range_tube (2,2) = 99999

      opt_eir%type         (3) = 2.0
      opt_eir%npts         (3) = -90000
      opt_eir%flux         (3) = 1.0
      opt_eir%flux_fraction(3) = 1.0
      opt_eir%species      (3) = 4
      opt_eir%species_index(3) = 1
      opt_eir%sorene       (3) = 0.0
      opt_eir%target       (3) = 1
      opt_eir%txtsou       (3) = 'default volume recombination'

      opt_eir%time  = 30
      opt_eir%niter = 0

      opt_eir%ntime = 0
      opt_eir%dtimv = 100.0E-06
      opt_eir%time0 = 0.0

      opt_eir%ilspt = 0

      opt_eir%tet_iliin = 2  ! Absorbing surface

      opt_eir%geom  = 2
      opt_eir%data  = 1
      opt_eir%trim  = 1
      opt_eir%alloc = 1.0

      opt_eir%photons = 0
      opt_eir%opacity = 0
      opt_eir%bgk     = 0

      opt_eir%ntorseg = 30
      opt_eir%torfrac = 1.0
      opt_eir%mat1    = 2
      opt_eir%mat2    = 2
      opt_eir%ctargt  = 300.0
      opt_eir%cwallt  = 300.0

      opt_eir%sur_n   = 0  ! Surface properties 
      opt_eir%tet_n   = 0  ! Tetrahedron mesh definition
      opt_eir%nadspc  = 0  ! Energy spectra definitions

      opt_eir%spcvx = 0.0  
      opt_eir%spcvy = 0.0
      opt_eir%spcvz = 0.0
      opt_eir%spc_p1 = -999.0
      opt_eir%spc_p2 = -999.0

      opt_eir%whipe   = 0  ! Debugging mode where the plasma density is set to very low values everywhere

      eirfp = 88     
      geofp = 88

c...  SOLPS related variables:
      solps_opt = 0
      nsolps_data = 0
      IF (ALLOCATED(solps_data)) DEALLOCATE(solps_data)
      IF (ALLOCATED(map_divimp)) DEALLOCATE(map_divimp)
      IF (ALLOCATED(solps_cen )) DEALLOCATE(solps_cen )
      IF (ALLOCATED(map_osm   )) DEALLOCATE(map_osm   )

      nhistory = 0

c...  Divimp options:
      opt_div%rib_n = 0
      opt_div%rib_format = -1
      opt_div%rib_r1 = -1.0E+6
      opt_div%rib_r2 =  1.0E+6
      opt_div%rib_z1 = -1.0E+6
      opt_div%rib_z2 =  1.0E+6
      opt_div%rib_rad_opt = 0
      opt_div%rib_pol_opt = 0

      opt_div%rib_pol_n_def = 11 
      opt_div%rib_pol_a_def = 0.3
      opt_div%rib_pol_b_def = 0.1
      opt_div%rib_pol_c_def = 1.0
      opt_div%rib_pol_d_def = 0.1

c...  User:
      CALL User_InitializeOptions

      RETURN
 99   STOP
      END
c
c ======================================================================
c
