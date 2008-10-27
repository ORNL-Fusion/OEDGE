c     -*-Fortran-*-
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

c        WRITE(0,*) '  BUFFER:',TRIM(buffer),mode

        IF (LEN_TRIM(buffer).EQ.0) THEN
c...      Comment or blank line, so continue:
        ELSEIF (buffer(1:1).EQ.'{'.OR.buffer(2:2).EQ.'{') THEN
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
      LOGICAL FUNCTION GetLine_Old(fp,buffer)
      IMPLICIT none

      INTEGER  , INTENT(IN)  :: fp
      CHARACTER, INTENT(OUT) :: buffer*(*)

      INTEGER i

      GetLine_Old = .TRUE. 

      DO WHILE(.TRUE.)
        READ(fp,'(A)',END=10) buffer

c        WRITE(0,*) 'BUFFER:',buffer(1:50)

c...    Remove portion of line after comment charcter:
        DO i = 1, LEN_TRIM(buffer)
          IF (buffer(i:i).EQ.'$'.OR.buffer(i:i).EQ.'*') EXIT
        ENDDO 
        buffer(i:LEN_TRIM(buffer)) = ' '

c...    Tag line discovered so backup file position and return
c       to main input file scanning loop:
        IF (buffer(2:2).EQ.'{') THEN
          GetLine_Old = .FALSE.
          BACKSPACE(fp)
          EXIT
        ELSEIF (LEN_TRIM(buffer).EQ.0) THEN
c...      Comment or blank line, so continue:
        ELSE
c...      Data line found:
          EXIT
        ENDIF
      ENDDO
c      WRITE(0,*) 'BUFFER:',buffer(1:50),GetLine_Old
      RETURN

 10   GetLine_Old = .FALSE.
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
      SUBROUTINE LoadOptions
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_legacy
      IMPLICIT none

      LOGICAL GetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

      INTEGER   fp,i,j
      CHARACTER buffer*1024,cdum1*512

      LOGICAL :: status = .TRUE.


c...  Set default values:
      CALL InitializeOptions


c...  Open log file (may be closed below if logfp is set to 0):
c      OPEN(UNIT=logfp,FILE='osm_log.dat',ACCESS='SEQUENTIAL',
c     .     STATUS='REPLACE',ERR=98)     


c...  Open OSM input file that contains option settings:
      fp = 99
      OPEN(UNIT=fp,FILE='osm.input',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)

c...  Scan input file looking for tags which are marked by a '{' in 
c     column 2 on an input line and contained by curly brackets:
      DO WHILE(GetLine(fp,buffer,WITH_TAG))

c...    Isolate tag string:
        DO i = 2, LEN_TRIM(buffer)
          IF (buffer(i:i).EQ.'}') EXIT
        ENDDO

        CALL ProcessInputTag(fp,i,buffer,status)

        IF (.NOT.status) EXIT
      ENDDO
 10   CONTINUE

      CLOSE(fp)

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

      INTEGER   fp,itag
      LOGICAL   status
      CHARACTER buffer*(*)

      WRITE(logfp,*) 'TAG:>'//buffer(3:itag-1)//'<'

c...  Check for reserved group markers which are indicated by the
c     first two characters in the tag, otherwise treat the tag as
c     miscellaneous:
      SELECTCASE (buffer(3:4))
        CASE ('C ')
          CALL LoadControlOption(fp,buffer,itag)
        CASE ('E ')
          CALL LoadEireneOption(fp,buffer,itag)
        CASE ('F ')
          CALL LoadFileOption(fp,buffer,itag)
        CASE ('S ')
          CALL LoadSourceOption(fp,buffer,itag)
        CASE ('T ')
          CALL LoadTransportOption(fp,buffer,itag)
        CASE DEFAULT
          CALL LoadMiscOption(fp,buffer,itag,status)
      ENDSELECT


      RETURN
 98   CALL ER('LoadOptions','Error reading file',*99)
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

      CHARACTER cdum1*512

      SELECTCASE (buffer(3:itag-1))
        CASE('C KIN')
          READ(buffer,*) cdum1,opt%nflukin
        CASE('C EIR')
          READ(buffer,*) cdum1,opt%eirene
        CASE('C EIR_TIME')
          READ(buffer,*) cdum1,opt_eir%time
      CASE DEFAULT
          CALL User_LoadOptions(fp,itag,buffer)
      ENDSELECT 

      RETURN
 99   STOP
      END

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

      LOGICAL   GetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

      INTEGER   i1
      REAL      stratum_type,version,rdum(7)

      SELECTCASE (buffer(3:itag-1))
        CASE('E NEUTRAL SOURCES')
          nstrata = 0
          DO WHILE(GetLine(fp,buffer,NO_TAG))
            nstrata = nstrata + 1
c            WRITE(0,*) 'BUFFER:',TRIM(buffer)
            READ(buffer,*) 
     .        strata(nstrata)%type,
     .        strata(nstrata)%npts,
     .        strata(nstrata)%flux,
     .        strata(nstrata)%flux_fraction,
     .        strata(nstrata)%species,
     .        strata(nstrata)%species_index,
     .        strata(nstrata)%sorene
            IF     (strata(nstrata)%type.EQ.1.0) THEN
              STOP 'NOT READY'
            ELSEIF (strata(nstrata)%type.EQ.2.0) THEN
              STOP 'NOT READY'
            ELSEIF (strata(nstrata)%type.EQ.3.0) THEN
              READ(buffer,*) rdum(1:7),
     .          strata(nstrata)%sorcos,
     .          strata(nstrata)%sormax,
     .          strata(nstrata)%sorad(1:6),
     .          strata(nstrata)%txtsou
            ELSE
              CALL ER('LoadEireneOption','Unknown stratum type',*99)
            ENDIF
          ENDDO            
c          WRITE(0,*) 'NSTRATA:',nstrata,rdum(1:6)
c          WRITE(0,*) 'NSTRATA:',nstrata,strata(nstrata)%sorad
c          STOP
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
      USE mod_sol28_global
      IMPLICIT none

      INTEGER  , INTENT(IN) :: fp,itag
      CHARACTER, INTENT(IN) :: buffer*(*)

      CHARACTER cdum1*512

      SELECTCASE (buffer(3:itag-1))
        CASE('F LOG')
          READ(buffer,*) cdum1,log
          IF (log.LT.0) THEN
            CLOSE(logfp)
            logfp = 0
            log   = -log
          ENDIF
          opt%log   = log
          opt%logfp = logfp
        CASE('F OSM_LOAD')
          opt%osm_load = 1
          WRITE(opt%f_osm_load,'(512X)')
          READ(buffer,*) cdum1,opt%f_osm_load
        CASE('F OSM_DIR')
          WRITE(opt%f_osm_dir,'(512X)')
          READ(buffer,*) cdum1,opt%f_osm_dir
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

      LOGICAL   GetLine
      INTEGER, PARAMETER :: WITH_TAG = 1, NO_TAG = 2

      INTEGER   i1
      LOGICAL   node_fit,node_data,ldum1
      REAL      node_type,rdum1
      CHARACTER cdum1
      TYPE(type_node) :: node_tmp

      status = .TRUE.

      SELECTCASE (buffer(3:itag-1))
        CASE('088')
          tarninter(HI) = 0
          DO WHILE(GetLine(fp,buffer,NO_TAG))
            tarninter(HI) = tarninter(HI) + 1
            READ(buffer,*) tarinter(tarninter(HI),1:4,HI)
          ENDDO
        CASE('089')
          tarninter(LO) = 0
          DO WHILE(GetLine(fp,buffer,NO_TAG))
            tarninter(LO) = tarninter(LO) + 1
            READ(buffer,*) tarinter(tarninter(LO),1:4,LO)
          ENDDO
        CASE('S74')
          READ(buffer(itag+2:itag+4),*) opt%s28mode  
          IF    (opt%s28mode.EQ.4.0) THEN
            osmns28 = 0            
            DO WHILE(GetLine(fp,buffer,NO_TAG))
              osmns28 = osmns28 + 1
              READ(buffer,*) osms28(osmns28,1:12)
            ENDDO
          ELSEIF (opt%s28mode.EQ.4.1) THEN
            node_fit = .FALSE.
            osmnnode = 0            
            DO WHILE(GetLine(fp,buffer,NO_TAG))
              CALL PadBufferR(15,buffer,0.0)
              osmnnode = osmnnode + 1
              READ(buffer,*) node_type
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
                  ldum1 = GetLine(fp,buffer,NO_TAG)
                  CALL PadBufferR(15,buffer,1.0)
                  READ(buffer,*) node_tmp%file_name,
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
c...            Load data for radial fits:
                osmnode(osmnnode)%type = 0.0
                READ(buffer,*) rdum1,
     .            osmnode(osmnnode)%fit_type,
     .            osmnode(osmnnode)%fit_psin(1:2),
     .            osmnode(osmnnode)%fit_shift,
     .            osmnode(osmnnode)%fit_quantity,
     .            osmnode(osmnnode)%fit_p(1:6)
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

c          DO i1 = 1, osmnnode
c            WRITE(0,*) 'nodes:',osmnode(i1)%type 
c          ENDDO

        CASE('030')
        CASE('E16')
        CASE('999')
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
      SUBROUTINE InitializeOptions
      USE mod_sol28_global
      USE mod_eirene06
      USE mod_legacy
      IMPLICIT none

      CALL InitializeLegacyVariables

c...  Output options:
      opt%log = 0
      log     = 1
      logfp   = 88

      opt%osm_load = 0

c...  Node interpolation data:
      opt%s28mode = 4.0
      osmnnode = 0

c...  

c...  OSM options:
      nion = 1

      opt%pin_data   = .FALSE.
      opt%p_ion      = 3.1
      opt%p_ion_frac = 100.0
      opt%p_rec      = 0.1 
      opt%p_ano      = 2
      opt%m_mom      = 0.1
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

c...  Eirene options:
      nstrata = 0

      opt_eir%time  = 30
      opt_eir%niter = 0

      opt_eir%dtimv = 0.0

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


      CALL User_InitializeOptions


      RETURN
 99   STOP
      END
c
c ======================================================================
c
