c     -*-Fortran-*-
c
c ======================================================================
c
      LOGICAL FUNCTION osmGetLine(fp,buffer,mode)
      USE mod_sol28_global
      IMPLICIT none

      INTEGER  , INTENT(IN)  :: fp,mode
      CHARACTER, INTENT(OUT) :: buffer*(*)

      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2, ALL_LINES = 3

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
      USE mod_sol28_global
      USE mod_legacy
      IMPLICIT none

      LOGICAL osmGetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

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
      USE mod_sol28_global
      IMPLICIT none

      INTEGER  , INTENT(IN) :: fp,itag
      CHARACTER, INTENT(IN) :: buffer*(*)

      LOGICAL osmGetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2, ALL_LINES = 3

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
      SUBROUTINE LoadEireneOption(fp,buffer,itag)
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_legacy
      USE mod_eirene06
      IMPLICIT none

      INTEGER  , INTENT(IN)  :: fp,itag
      CHARACTER  :: buffer*(*)

      LOGICAL   osmGetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

      INTEGER   i1,idum(5)
      CHARACTER cdum*1024
      REAL      stratum_type,version,rdum(7)

      SELECTCASE (buffer(3:itag-1))
c       ----------------------------------------------------------------
        CASE('EIR IMPURITY SPUTTERING')
          CALL ReadOptionI(buffer,1,opt_eir%ilspt) 
c       ----------------------------------------------------------------
        CASE('EIR VOID GRID')
          opt_eir%nvoid = 0
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
            opt_eir%nvoid = opt_eir%nvoid + 1
c            WRITE(0,*) 'BUFFER:',TRIM(buffer)
            READ(buffer,*) idum(1)
            IF (idum(1).EQ.-2) THEN
              READ(buffer,*) 
     .          opt_eir%void_zone(    opt_eir%nvoid),
     .          opt_eir%void_grid(1:2,opt_eir%nvoid)
            ELSE
              READ(buffer,*) 
     .          opt_eir%void_zone(    opt_eir%nvoid),
     .          opt_eir%void_grid(1:2,opt_eir%nvoid),
     .          opt_eir%void_wall(1:2,opt_eir%nvoid),
     .          opt_eir%void_add (1:2,opt_eir%nvoid),
     .          opt_eir%void_res (    opt_eir%nvoid),
     .          opt_eir%void_hole(1:2,opt_eir%nvoid),
     .          opt_eir%void_code(    opt_eir%nvoid),
     .          opt_eir%void_ne  (    opt_eir%nvoid),
     .          opt_eir%void_te  (    opt_eir%nvoid),
     .          opt_eir%void_ti  (    opt_eir%nvoid)
            ENDIF
          ENDDO
c       ----------------------------------------------------------------
        CASE('E NEUTRAL SOURCES')
          opt_eir%nstrata = 0
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
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
     .              opt_eir%type(opt_eir%nstrata).EQ.3.1) THEN  ! Point source injection (gas puff, beams)
              READ(buffer,*) rdum(1:7),
     .          opt_eir%sorcos   (opt_eir%nstrata),
     .          opt_eir%sormax   (opt_eir%nstrata),
     .          opt_eir%sorad(1:6,opt_eir%nstrata),
     .          opt_eir%txtsou   (opt_eir%nstrata)
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
     .        opt_eir%ispsrf    (opt_eir%nadspc),  ! Surface index, <0=non-default standard, >0=additional surfaces          
     .        opt_eir%ispsrf_ref(opt_eir%nadspc),  ! Which code does the surface index refer to?                             
     .        opt_eir%iptyp     (opt_eir%nadspc),  ! Species type eg 1=atoms, 2=molecules, 3=test ions, 4=?                  
     .        opt_eir%ipsp      (opt_eir%nadspc),  ! Species sub-index eg, 1=first atom species, 2=second atom species, etc.   
     .        opt_eir%isptyp    (opt_eir%nadspc),  ! Spectrum type wrt units, 1=1/eV/s, 2=1/s                                  
     .        opt_eir%nsps      (opt_eir%nadspc),  ! Number of bins                                                          
     .        opt_eir%spcmn     (opt_eir%nadspc),  ! Lower bound of energy range for spectrum                                
     .        opt_eir%spcmx     (opt_eir%nadspc),  ! Upper bound                                                             
     .        opt_eir%idirec    (opt_eir%nadspc)   ! If >0 then a projection on a direction is used in the statistics (??)   
            IF (opt_eir%idirec(opt_eir%nadspc).NE.0) 
     .        READ(buffer,*) idum(1:2),cdum   ,idum(1:4),
     .                       rdum(1:2),idum(1),
     .          opt_eir%spcvx(opt_eir%nadspc),     ! Don't know really, but it was in the example that AK sent, originally from VK
     .          opt_eir%spcvy(opt_eir%nadspc),
     .          opt_eir%spcvz(opt_eir%nadspc)
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
      USE mod_options
      USE mod_legacy
      USE mod_eirene06
      IMPLICIT none

      INTEGER  , INTENT(IN)  :: fp,itag
      CHARACTER  :: buffer*(*)

      LOGICAL   osmGetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

      INTEGER   i1
      REAL      stratum_type,version,rdum(7),length(2)

      SELECTCASE (buffer(3:itag-1))
        CASE('FIL OPT')
          CALL ReadOptionI(buffer,1,opt_fil%opt)
        CASE('FIL TARGET FLUX')
          CALL ReadOptionI(buffer,1,opt_fil%target_flux)
        CASE('FIL LENGTH')
          WRITE(0,*) 'BUFFER:',TRIM(buffer)
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
          WRITE(0,*) 'opt%f_eirene_15:',TRIM(opt_eir%f_eirene_15)
        CASE('GRID FORMAT')
          CALL ReadOptionI(buffer,1,opt%f_grid_format)
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
      USE mod_sol28_global
      IMPLICIT none

      INTEGER  , INTENT(IN) :: fp,itag
      CHARACTER, INTENT(IN) :: buffer*(*)

      LOGICAL osmGetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

      SELECTCASE (buffer(3:itag-1))
        CASE('SOL APPLICATION')
          opt%sol_n = 0
          DO WHILE(osmGetLine(fp,buffer,NO_TAG))
            opt%sol_n = opt%sol_n + 1
            READ(buffer,*) opt%sol_tube  (1:2,opt%sol_n),
     .                     opt%sol_option(    opt%sol_n)
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
      SUBROUTINE LoadSOLPSOption(fp,buffer,itag)
      USE mod_solps_params
      USE mod_solps
      IMPLICIT none

      INTEGER  , INTENT(IN) :: fp,itag
      CHARACTER, INTENT(IN) :: buffer*(*)

      CHARACTER cdum1*1024

      LOGICAL osmGetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

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
      USE mod_sol28_wall
      IMPLICIT none

      INTEGER  , INTENT(IN) :: fp,itag
      CHARACTER, INTENT(IN) :: buffer*(*)

      CHARACTER cdum1*1024

      LOGICAL osmGetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

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
      USE mod_sol28_global
      USE mod_legacy
      IMPLICIT none

      INTEGER  , INTENT(IN)  :: fp,itag
      CHARACTER  :: buffer*(*)
c      CHARACTER, INTENT(IN)  :: buffer*(*)
      LOGICAL  , INTENT(OUT) :: status

      LOGICAL osmGetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

      INTEGER   i1
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
     .            osmnode(osmnnode)%potential
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
     .                  osmnode(osmnnode)%fit_quantity,
     .                  osmnode(osmnnode)%fit_p(1:9)
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
      USE mod_legacy
      USE mod_solps
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
      opt%p_ion      = 3
      opt%p_ion_exp  = 0.1
      opt%p_ion_frac = 100.0
      opt%p_rec      = 0 
      opt%p_ano      = 2
      opt%m_mom      = 0
      opt%m_fit      = 2
      opt%m_ano      = 2
      opt%m_ano_dist = 1
      opt%m_ano_exp  = 0.0

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
      opt_eir%nstrata = 0

c      opt_eir%nvoid = 0
      opt_eir%nvoid = 1
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

      opt_eir%time  = 30
      opt_eir%niter = 0

      opt_eir%ntime = 0
      opt_eir%dtimv = 100.0E-06
      opt_eir%time0 = 0.0

      opt_eir%ilspt = 0

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

      opt_eir%nadspc  = 0  ! Energy spectra definitions

      eirfp = 88     
      geofp = 88

c...  SOLPS related variables:
      solps_opt = 0
      nsolps_data = 0
      IF (ALLOCATED(solps_data)) DEALLOCATE(solps_data)
      IF (ALLOCATED(map_divimp)) DEALLOCATE(map_divimp)
      IF (ALLOCATED(solps_cen )) DEALLOCATE(solps_cen )
      IF (ALLOCATED(map_osm   )) DEALLOCATE(map_osm   )

c...  User:
      CALL User_InitializeOptions

      RETURN
 99   STOP
      END
c
c ======================================================================
c
