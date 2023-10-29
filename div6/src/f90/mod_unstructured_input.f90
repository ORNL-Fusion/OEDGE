module unstructured_input


  !     This module contains code related to the unstructured input options 

  !
  ! ======================================================================
  !
  ! contains all the initialization and read code for DIVIMP unstructured input except the SOL22 block which is grouped
  ! with the SOL22 code for standalone operation and use in LIM. 

  implicit none


    integer,parameter,private :: maxlen=1024

    integer, PARAMETER :: MAXTAG = 1000
    integer :: ntaglist = 0
    character*3 :: taglist(maxtag) 
    ! jdemod - inputflag seems to be something old that isn't in use anymore
    !INTEGER     :: inputflag(100)

    logical :: read_nimbin = .false.
    
    ! jdemod - this is for old code that appears to be permanently deactivated - but used in two routines here - so move to module
    CHARACTER*64     machine2
  
    integer :: nqs  ! Local variable used to read in dwell times - this must be equal to nizs+1 for a time dependent case
    ! jdemod - these are for reading in debugging steps
    integer :: istepl  ! local for loading cstepl (real) using RDI (?) 
    integer :: istepn  ! local for loading cstepn (real) using RDI (?) 
    integer :: istepv  ! local for loading cstepv (real) using RDI (?) 
    
    
contains


  

  subroutine scan_input_file
    use mod_params
    use mod_io_units
    use mod_reader
    use mod_io
    use mod_dynam4
    implicit none

    ! Update some parameters based on input file contents
    ! Also attempt to determine whether the input file is structured or tagged

    integer :: nizs_tmp, cion_tmp, nimps_tmp, nimps2_tmp, nts_tmp, nitersol_tmp, nfla_tmp
    real,allocatable :: dwelfs_tmp(:)
    !logical :: is_structured

    character*(maxlen) :: line
    character*3 :: tag
    character*2 :: tag_tn
    character*1 :: tag_mod
    integer :: tn_count

    integer :: ierr
    logical :: err

    ! jdemod
    ! - access the input file
    ! - read in any updates to the parameter values from their default values
    ! - determine the type of input file - tagged or untagged 

    ! parameter tags start with the # sign and are in series Z

    ! Important tags for determining paramters:
    !
    ! 1) MAXIZS (No Z tag)
    ! Ideally the following are the same. If not then take the larger. This defines parameter MAXIZS if different from the default. 
    ! If only one is specified then the maximum between the input and default is chosen. The default values for both are "6" for Carbon impurity.
    !
    ! I15 Maximum ionization state to be followed
    ! S05 Atomic number of impurity ions
    !
    ! 2) MAXIMP (No Z tag)
    !
    ! MAXIMP is the sum of S11 and S12
    !
    ! S11 Number of impurity ions to follow
    ! S12 Number of supplemental impurity ions to follow
    !
    ! 3) MAXNTS (minimum of NTS or 1)
    !
    ! S17 Dwell Time Factors for time-dependent analysis
    !
    ! 4) MAXPINITER (Maximum number of SOL/PIN iterations)
    !
    ! P39 Number of PIN/SOL iterations
    !
    ! 5) MAXNFLA (Max number of fluid species in fluid code data file)
    !
    ! F05 Number of fluids in B2 background
    !
    ! Needs to be a minimum of 1 due to nfla2 used in eirediv.f which is a local parameter=1
    !

    
    ierr = 0
    is_structured = .false.
    read_nimbin = .false. 
    
    tn_count = 0

    nizs_tmp = 0
    cion_tmp = 0
    nimps_tmp = 0
    nimps2_tmp = 0
    nts_tmp = 0
    nitersol_tmp = 0
    nfla_tmp = 0
    
    err = .false.
    
    do while (.not.err)

       ! loop through input file

       call divrd(ierr,line)
      
       if (ierr.eq.0) then
          ! check for parameter tags - series Z
          
          WRITE(TAG,'(A3)') LINE(3:5) 
          tag_mod = line(2:2)
          tag_tn  = line(2:3)
          ! jdemod - add support for case insensitive tags
          call upcase(tag)
          call upcase(tag_tn)
          
          ! possible parameter over ride found - parameters use series Z tags
          if (tag_mod == '#') then 
          
             if (tag(1:3).eq.'Z01') then 
                ! Z01 : MAXNRS - Maximum number of rings on the grid
                CALL divrd(maxnrs,.TRUE.,100,.false.,0,'Maximum number of rings on the grid',IERR)

             elseif (tag(1:3).eq.'Z02') then 
                ! Z02 : MAXNKS - Maximum number of knots on each ring of the grid
                CALL divrd(maxnks,.TRUE.,200,.false.,0,'Maximum number of knots on each ring of the grid',IERR)

             elseif (tag(1:3).eq.'Z03') then 
                ! Z03 : MAXPTS - max number of wall data/SOL23 info/ADAS variables/etc
                CALL divrd(maxpts,.TRUE.,500,.false.,0,'max number of wall data/SOL23 info/ADAS variables/etc',IERR)

             elseif (tag(1:3).eq.'Z04') then 
                ! Z04 : MAXSEG - Max number of wall segments (used for wall definition and wall fluxes)
                CALL divrd(maxseg,.TRUE.,500,.false.,0,'Max number of wall segments (used for wall definition and wall fluxes)',IERR)                

             elseif (tag(1:3).eq.'Z05') then 
                ! Z05 : MAXNWS - number of walk steps recorded
                CALL divrd(maxnws,.TRUE.,500,.false.,0,'Number of walk steps recorded',IERR)
                
             elseif (tag(1:3).eq.'Z06') then                 
                ! Z06 : MSOLPT - Number of points for detailed SOL backgrounds - soledge
                CALL divrd(msolpt,.TRUE.,100,.false.,0,'Number of points for detailed SOL backgrounds - soledge',IERR)
                
             elseif (tag(1:3).eq.'Z07') then 
                ! Z07 : MAXDATX - Max number of experimental data points allowed to be loaded
                CALL divrd(maxdatx,.TRUE.,500,.false.,0,'Max number of experimental data points allowed to be loaded',IERR)
                
             endif
          elseif (tag_mod == '+') then 
             ! Look for specific tags needed to set parameters correctly
             if (tag(1:3).eq.'I15') then 
                ! I15 : NIZS - Maximum ionization state to be followed
                CALL divrd(NIZS_tmp,  .TRUE.,  0, .FALSE.,MAXIZS,'MAX IZ STATE',   IERR)

             elseif (tag(1:3).eq.'S05') then 
                ! S05 : CION - Atomic number of impurity ions
                call divrd(CION_tmp,  .TRUE. , 1 ,.FALSE., 0 ,'IMP ATOMIC NUMBER', IERR)

             elseif (tag(1:3).eq.'S11') then 
                ! S11 : NIMPS - Number of particles to follow
                call divrd(NIMPS_tmp, .TRUE.,  1, .FALSE.,0,'NO OF IONS',     IERR)

             elseif (tag(1:3).eq.'S12') then 
                ! S12 : NIMPS2 - Number of supplemental particles to follow
                call divrd(NIMPS2_tmp,.TRUE.,0,.FALSE.,0,'NUM SUP IONS',IERR)

             elseif (tag(1:3).eq.'S17') then 
                ! S17 : NTS - number of walk steps recorded
                call divrda(DWELFS_tmp,NTS_tmp,MAXNTS,  0.0,machhi,.TRUE.,'DWELL Time FACTORS',IERR)
                
             elseif (tag(1:3).eq.'P39') then                 
                ! P39 : NITERSOL - Number of OSM iterations with PIN (EIRENE or NIMBUS)
                call divrd(NITERSOL_tmp,.TRUE.,0 ,.FALSE.,maxpiniter,'NO. OF PIN ITER. ',IERR)
                
             elseif (tag(1:3).eq.'F05') then 
                ! F05 : NFLA - Number of fluids in B2 solution file. 
                call divrd(nfla_tmp,.TRUE.,1,.false.,maxnfla,'NUM.FLUIDS IN B2 BG',IERR)
                
             elseif (tag(1:3).eq.'A01') then 
                ! check series of tags that are likely essential to most cases - if any of these are present then
                ! it is likely a tagged input file.
                ! A01 - title
                ! A02 - run comment 
                ! S05 - atomic number of impurity
                ! S11 - number of impurities to follow
                
             endif

          elseif (tag_tn == 'TN') then
             ! DIVIMP input files have a lot of lines containing a TN designation referring to a TNote. Untagged files
             ! will have many lines that start with TN. If there are more than 10 (arbitrary) then the file is a structured
             ! and untagged input file. 
             tn_count = tn_count +1
          elseif (index(line,'&NIMBIN') .ne. 0) then
             read_nimbin = .true.
          endif
             
       else
          err = .true.
       endif
          
    end do

    ! rewind input file
    call divrd

    ! Is the input file structured or unstructured - based on the number of TN references in column 2:3 of the input line.
    if (tn_count > 10) then 
       is_structured = .true.
       write(0,*) 'SCAN_INPUT_FILE: READING STRUCTURED INPUT FILE'
    else
       ! is_structured defaults to false so setting it here isn't strictly necessary
       is_structured = .false.
       write(0,*) 'SCAN_INPUT_FILE: READING TAGGED INPUT FILE'       
    endif

    
    !     jdemod - After NIZS and CION have been read in - a value can be
    !              assigned to maxizs and used for storate allocation      
    !
    !     if nizs is not a tagged input then leave maxizs at the default value
    !     The default value for CION is 6 for carbon - use 6 as a minimum for maxizs even if nizs specified smaller.
    !     The only cost is unused storage. 
    !  
    !     If CION is specified then this takes precedence since it is the atomic number of the impurity and represents
    !     the maximum number of charge states to be followed. 
    !     

    ! MAXIZS
    
    if (cion_tmp > 0) then
       maxizs   = cion_tmp
       if (nizs_tmp > cion_tmp) then
          call errmsg('UNSTRUCTURED_INPUT:SCAN_INPUT_FILE: SPECIFIED NIZS (TAG I15) IS GREATER THAN IMPURITY SPECIES ATOMIC NUMBER (TAG S05):',cion_tmp)
          stop 'Invalid NIZS: TAG I15'
       endif
    else       
       ! NOTE: It is an error if NIZS is greater than CION
       ! The default value if CION is 6 for carbon if cion_tmp not specified
       if (nizs_tmp > 6) then
          call errmsg('UNSTRUCTURED_INPUT:SCAN_INPUT_FILE: SPECIFIED NIZS (TAG I15) IS GREATER THAN IMPURITY SPECIES ATOMIC NUMBER DEFAULT (TAG S05):',6)
          stop 'Invalid NIZS: TAG I15'
       endif
    endif

    ! MAXIMP
    
    if ((nimps_tmp + nimps2_tmp) > 0 ) then 
       ! only over ride maximps if input values are specified
       ! NOTE: MAXIMP is recalculated at the start of allocate_storage_div using the same formula - this is necessary for
       ! untagged cases where these values are not known until after the input file is read in. 
       maximp = nimps_tmp+nimps2_tmp+1
    endif

    ! MAXNTS
    
    if (nts_tmp > 0) then
       maxnts = max(nts_tmp,1)
    endif

    ! MAXPINITER
    
    if (nitersol_tmp > 0) then
       maxpiniter = nitersol_tmp
    endif
    
    ! MAXNFLA

    if (nfla_tmp > 0) then
       maxnfla = nfla_tmp
    endif

    ! - this can be put back into regular input file allocation now that maxizs is derived from the input file
    call allocate_mod_dynam4_input_special ! check if this is required 


    !write(0,*) 'SCAN_INPUT_FILE:',is_structured,tn_count,maxizs,maxnfla,maxpiniter,maximp,maxnts,read_nimbin

    
  end subroutine scan_input_file


  subroutine read_input_file(ierr)
    use mod_params
    use mod_io
    use mod_pindata
    implicit none
    integer :: ierr
    logical :: err
    character*(maxlen) :: line
    
    ! This routine is the driver for reading an unstructured input file - only the tags matter and are processed by the ReadUnstructuredInput routine
    ! The code here cycles through the input file exiting when all lines have been read. 

    ! Any NIMBUS input is expected only at the end of the input file. All lines from the one including &NIMBIN to the end of the file are read into the
    ! NIMBUS input block. This should only be included IF NIMBUS is used (which is pretty much never these days)
    
    err = .false.
    
    do while (.not.err)

       ! loop through input file

       call divrd(ierr,line)

       if (line(1:1).eq.'$'.or.line(1:1).eq.'c'.or.line(1:1).eq.'C'.or.(trim(line)==''.and.ierr==0)) cycle   ! ignore comment lines in the input file - also ignore blank lines

       if (ierr.eq.0) then
          if (read_nimbin.and.index(line,'&NIMBIN').ne.0) then 
          
             ! NIMBUS input namelist - not used - from &NIMBIN to &END if found in the input file
             ! NIMBUS reads the input file directory to get this info - so reading it in DIVIMP is
             ! purely for documentation purposes
             !+H14 Nimbus input namelist: (see, for example, PINPGX for definitions)
             CALL RDCAR(CNIMBIN,nlines,maxlines,'NIMBUS NAMELIST INPUT',IERR)
             ! This is the last part of the input file read. If the NIMBUS namelist is present it is at the end
             ! After this returns set err = .true. and ierr = 0
             ! This is necessary since reading another line of input at this point is actually trying to read
             ! beyond the EOF which has already been encountered. Rather than EOF this returns a nonstandard error code.

             err = .true.
             ierr = 0
          else
             call readunstructuredinput(line)
          endif
       else
          err = .true.
          ! if the error is end of file then reset ierr to zero since this is the expected normal exit
          ! the fortran 2003 standard only specifies that EOR and EOF must be less than zero. 
          if (ierr.lt.0) then 
             ierr = 0
          endif
       endif
          
    end do

  end subroutine read_input_file

  
  subroutine after_read_input
    use mod_params
    use mod_comtor
    use mod_diagvel
    use mod_pindata
    use mod_dynam4
    use mod_cedge2d
    use ero_interface
    use mod_fperiph_com
    use comhc
    use mod_rundiv_local
    use mod_slcom
    !use mod_div_input
    use mod_lambda
    use mod_promptdep
    use mod_sol22_input
    implicit none

    integer :: ierr, in
    real :: lpdat_conv
    real :: zo
    
    call ValidateUnstructuredInput

    call verify_sol22    

       !     If using this model, then we must ensure CENUTD = 0, it doesn't 
       !     make sense otherwise. 
       if (cneutd /= 0.and.csputopt == 8) then 
          write(0,*) 'Warning! SiC model in use, CNEUTD forced to 0.' 
          write(0,*) 'Change CNEUTD to 0 to avoid this warning.' 
          cneutd = 0 
       end if

       
    
    ! P02
    ! move to after input file read
    if (ccoreopt.gt.6.and.ccoreopt.ne.28) then
       write(0,*) 'ERROR READIN: INVALID CORE OPTION'
       ierr=1
       return
    endif

    ! P03

    ! move to after input file read

    ! jdemod - code needs to move to start of bgplasma likely - or at least after input file read
    !...  Need to store these because they are overwritten in BGPLASMA but are needed in GETMODEL:
    orgcioptf = cioptf
    orgcioptg = cioptg
    IF (CSECSOL.EQ.-2) CSECSOL = CIOPTF

    !     IF CNEUTH OR CNEUTI ARE -1 THEY NEED TO BE SET TO THE VALUES
    !     READ IN FOR CNEUTB AND CNEUTC.
    !
    ! move to AFTER the input file has been completely read in
    IF (NVAOPT.EQ.-1) NVAOPT = CNEUTC
    IF (CNEUTH.EQ.-1) CNEUTH = CNEUTB
    IF (CNEUTI.EQ.-1) CNEUTI = CNEUTC


    ! move to after input file read in
    if (cneutd2.eq.-1) cneutd2 = cneutd


    ! check this ...
    ! move to after input file read
    IF (CIOPTB.EQ.0.AND.(CIOPTC.EQ.0.or.cioptc.eq.4).AND.(CIOPTD.EQ.0.OR.CIOPTD.EQ.3)) IRSPEC=2*MAXNRS

    rizb = real(cizb)
    
    !     move to after input file read 
    !      
    !     Adjust the Target conditions - by multiplication factors
    !
    !     jdemod 
    !
    !     If the lpdat switch is set to 2 convert particles/s to Amperes for
    !     compatibility elsewhere in the code. 
    !      
    if (lpdatsw.eq.2) then 
       !
       !        Note: after conversion the data is then in A and the lpdatsw needs to be
       !        set accordingly. 
       !
       write(6,'(a)') 'INPUT: LPDATSW : Initial value 2 : Isats converted from particles/s to Amp : LPDATSW set to 1'
       !
       lpdat_conv = ech
       lpdatsw = 1
    else
       lpdat_conv = 1.0
    endif

    if (nlpdati.gt.0) then 

       do in = 1,nlpdati
          lpdati(in,2) =  lpdati(in,2) * te_mult_i
          lpdati(in,3) =  lpdati(in,3) * ti_mult_i
          lpdati(in,4) =  lpdati(in,4) *  n_mult_i * lpdat_conv
       end do
    endif

    if (nlpdato.gt.0) then 
       do in = 1,nlpdato
          lpdato(in,2) =  lpdato(in,2) * te_mult_o
          lpdato(in,3) =  lpdato(in,3) * ti_mult_o
          lpdato(in,4) =  lpdato(in,4) *  n_mult_o * lpdat_conv
       end do
    endif


    ! Setting region diffusion coefficients
    ! IF values < 0.0 then set to base Dperp

    ! PFZ/TRAPPED plasma
    if (cdperpt.le.0.0) cdperpt = cdperp

    ! Far Periphery
    IF (CDPERPFP.Le.0.0) CDPERPFP = CDPERP

    ! Core
    if (cdperpc.le.0) cdperpc = cdperp

    
    ! rescale normal angle for launches from degrees to radians
    CSNORM = CSNORM / RADDEG


    !
    !     RESET NIMPS2 - if an Ion injection has been specified and
    !                    not a neutral launch - set nimps2 = 0
    !
    ! move to after input file read
    if (cneuta.eq.1) nimps2 = 0


    CSTEPN = REAL (ISTEPN)
    CSTEPL = REAL (ISTEPL)
    CSTEPV = REAL (ISTEPV)

    ! jdemod - move setting of debug velocity flag to point where it
    ! is read in
    if (cstepv.ne.0.0) then
       debugv = .true.
    else
       debugv = .false.
    endif




    ! move to after input file read
    !
    !     PIN expects a value of -1 to generate a new seed while the DIVIMP
    !     documentation indicates <0 = 1, 0=generate new ... so we need
    !     to change the value of piniseed here to match the EIRENE
    !     expectations.(Other values are passed along as is). 
    !     
    if (piniseed.lt.0) then 
       piniseed = 1
    elseif (piniseed.eq.0) then 
       piniseed = -1
    endif

    ! move to after input file read - checks
    IF (NITERS.GT.1 .AND. (CIOPTB.EQ.3.or.cioptb.eq.9)) then
       CALL PRC ('READIN: ITERATIONS (TAG S20) > 1 INCOMPATIBLE WITH COLL OPT 3 OR 9 (TAG T02)')
       IERR = 1
       stop 'READIN: ITERATIONS (TAG S20) > 1 INCOMPATIBLE WITH COLL OPT 3 OR 9 (TAG T02)'
    endif

    ! ---------------   move to input read?
    IF (IMODE.EQ.1.AND.NQS.LT.NIZS+1) then
       CALL PRC ('READIN: NOT ENOUGH DWELL TIMES GIVEN')
       IERR = 1
       call errmsg('INPUT FILE ERROR: INSUFFICIENT DWELL TIMES (S16) GIVEN FOR TIME DEPENDENT CASE: NUMBER = NIZS+1 REQUIRED NIZS=',nizs)
       stop 'INSUFFICIENT DWELL TIMES: TAG S16'

       ! RETURN     ! need to decide appropriate response? Perhaps move this error message to point where input is read
    elseif (nqs.gt.0) then 
       DWELTS(-1) = DWELTS(0)
    endif


    ! move to after input file read
    ZO = REAL(CZO)
    CHZO = 1.56*(1.+1.41*ZO)*(1.+0.52*ZO)/((1.+2.65*ZO)*(1.+0.285*ZO))


    ! If option C01 is not set then consider establishing a default set of S values for the ion leakage diagnostic
    !
    ! TAG:  +C01    Initialization
    !
    ! Comments from input
    !
    ! Example
    !'TN982     Number of S-values  :-'  4
    !   5.0
    !   10.0
    !   15.0
    !   20.0

    !
    !     Set up the flag indicating that leak checking code should be active
    !
    if (cleaksn.gt.0) then
       checkleak = .true.
    else
       checkleak = .false.
    endif

    !
    !--------------------------------------------------------------------------
    ! 
    !     CHECKS AND VERIFICATION OF SOME INPUT VALUES:
    !
    !--------------------------------------------------------------------------   
    !
    !
    !     SOME INPUT DATA VERIFICATION
    !
    !
    !     If the target option specified requires the platco array to
    !     store data (as in target option 1), then the geometry data
    !     can not be loaded on top.
    !
    !    Note: if the target option is 3 or the wall option 4 or 5
    !           and cgeoopt has not been specified then an error
    !           has occurred.
    !
    IF (CGEOOPT.EQ.-1.AND.(ctargopt.eq.3.or.cneur.eq.3)) then
       call prc('Invalid Input ... target or wall option that')
       call prc('requires pre-loaded data has been specified ')
       call prc('but specific pre-loaded grid data has not   ')
       call prc('selected.')
       stop
    endif
    !
    !     Check to make sure that Trap option 3 and a double-null ITER
    !     grid are not specified at the same time - these options share
    !     a common array with different meanings at this time and are thus
    !     incompatible.
    !
    if (cgridopt.eq.2.and.ctrap.eq.3) then
       call prc('ITER grid is incompatible with TRAP option 3.')
       call prc('Please respcify input - program ending ...')
       stop
    endif
    !
    !     Error - make sure PIN is running if NIMBUS wall specified.
    !
    if (cpinopt.eq.0.and.(cneur.eq.5.or.ctrap.eq.5)) then
       !
       call prc('ERROR: Incompatible NEUTRAL Wall/Trap Wall Options')
       call prc('Neutral Wall Option 5 or Trap Wall Option 5')
       call prc('specified without corresponding call to PIN/NIMBUS')
       call prc('Program stopping')
       !
       !        Write to screen for easy notification
       !
       write(0,*) 'ERROR: Incompatible NEUTRAL Wall/Trap Wall Options'
       write(0,*) 'Neutral Wall Option 5 or Trap Wall Option 5'
       write(0,*) 'specified without corresponding call to PIN/NIMBUS'
       write(0,*) 'Program stopping'
       stop
    endif
    !
    !     Warning - if one of the SOL - Plasma Decay or Temperature Gradient
    !               options is specified as 99 and at least one is NOT
    !             - issue a Warning!
    !
    if (cioptg.eq.99.or.cioptf.eq.99.or.cioptk.eq.99.or.cioptl.eq.99) then

       if (cioptg.ne.99.or.cioptf.ne.99.or.cioptk.ne.99.or.cioptl.ne.99) then
          !
          !           Issue Warning - since mixed file data specified.
          !
          call prb
          call prc(' WARNING --- WARNING --- WARNING !!! ')
          call prb
          call prc(' POSSIBLY INCONSISTENT INPUT: ')
          call pri(' SOL OPTION          : ',cioptf)
          call pri(' PLASMA DECAY OPTION : ',cioptg)
          call pri(' Te GRADIENT OPTION  : ',cioptk)
          call pri(' Ti GRADIENT OPTION  : ',cioptl)
          call prb
          call prc(' AT LEAST ONE HAS BEEN SPECIFIED AS FILE INPUT (99)')
          call prc(' BUT NOT ALL - THIS MAY PRODUCE UNEXPECTED RESULTS!')
          call prb
       endif
    endif
    !
    !     e2dtargopt=2 requires flux data for each flux tube to be specified.
    !     e2dtargopt=3,4
    !
    if ((e2dtargopt.eq.2.or.e2dtargopt.eq.3.or.e2dtargopt.eq.4.or.e2dtargopt.eq.5)&
         .and.(fluxpts.le.0.and.readaux.ne.1)) then
       call prc('ERROR: Incompatible Input.')
       call prc('- Fluid code target condition option 2, 3 or 4 has been')
       call prc('  specified without corresponding fluxes for')
       call prc('  each target element - FLUXINFO - or using')
       call prc('  FC Target Option and Read Auxiliary Input to ')
       call prc('  extract down fluxes from an external file.')
       call prc('Program stopping')
       stop
    elseif (e2dtargopt.ne.0.and.cre2d.eq.0) then
       call prc('ERROR: Incompatible Input.')
       call prc('- A Fluid Code  target condition option has been')
       call prc('  specified without the corresponding option')
       call prc('  to also read the Fluid Code background for reference.')
       call prc('Program stopping')
       stop

    endif
    !
    !     Error - CIOPTI=8 is invalid except for Carbon impurity in a
    !             deuterium plasma.
    !
    if ((ciopti.eq.8.or.ciopti.eq.9).and.(crmb.ne.2.0.or.crmi.ne.12.0)) then
       call errmsg('INPUT ERROR:','Incompatible Input: CX Recombination option 8 is ONLY compatible with a carbon impurity in a deuterium plasma.')
       stop 'Invalid CX option'
    endif
    !
    !      IF Pinch option 4 is specified then a Probability Distribution
    !      Function must be loaded - if not - the pinch option is turned 
    !      OFF. 
    !
    if (pinchopt.eq.4.and.pinch_npdf.le.0) then 
       !
       call errmsg('INPUT ERROR:','PINCH OPTION 4 SELECTED BUT PDF NOT SPECIFIED IN INPUT :'//&
            ' PINCH OPTION TURNED OFF (SET TO 0)')
       pinchopt = 0
    endif

    !     Issue an error message if fp_neut_opt is non-zero for
    !     any periphery option except 3.        
    !     
    if (fp_neut_opt.ne.0.and.fpopt.ne.3) then
       write(0,'(a,i8,a,i8)')  'WARNING (ERROR): PERIPHERY IONIZATION OPTION ',&
            fp_neut_opt,' HAS BEEN SPECIFIED WITH AN INCOMPATIBLE PERIPHERY OPTION',fpopt
       write(6,'(a,i8,a,i8)') 'WARNING (ERROR): PERIPHERY IONIZATION OPTION ',&
            fp_neut_opt,' HAS BEEN SPECIFIED WITH AN INCOMPATIBLE PERIPHERY OPTION',fpopt
    endif


    !       
    !      CSTMAX is set in the rundiv main program
    !
    !      CSTMAX = 10.0 / QTIM
    !
    !
    !
    if (neut2d_opt.eq.2.and.ero_particle_launch_opt.eq.0) then 
       call prc('ERROR: Incompatible Input.')
       call prc('       Neut2d option = 2 specified without')
       call prc('       ERO particle launch being selected')
       call prc('       K40 - ERO particle launch option')
       call prc('       must be set and ERO particle data')
       call prc('       file must be available')
       neut2d_opt = 0

       write(0,'(a)') 'WARNING: NEUT2D option set to ERO'//&
            ' but ERO particle launch option (K40) turned'//&
            ' off. NEUT2D OPT set to 0'

    endif

    ! jdemod - hc_lambda_calc input option has been deprecated in
    ! favor of a global defintion of the coulomb logarithm in
    ! mod_lambda.f90
    ! 
    hc_lambda_calc = lambda_opt


      ! Prompt redeposition option 3 is an ERO-based scaling for W only.
      if (prompt_depopt.eq.3.and.cion.ne.74) then
         call errmsg('INPUT FILE ERROR: Prompt redeposition option 3 (TAG I07) only'//&
              'applicable to W. Change ion to W (TAG S05) or change prompt redeposition option.',&
              'PROGRAM STOPPING')
        stop 'INPUT ERROR: Prompt dep opt 3 (I07) only applicable to W (S05)'
      endif
    
    call pr_trace('READ_INPUT_FILE','END')

    RETURN
  end subroutine after_read_input


  SUBROUTINE InitializeVariables
    USE mod_divimp
    use mod_params
    use mod_slcom
    use mod_cgeom
    use mod_cadas
    use mod_comtor
    use mod_pindata
    use mod_cedge2d
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    use mod_reiser_com
    use mod_driftvel
    use mod_fperiph_com
    use mod_dperpz
    use mod_slcom_sol28
    use mod_sl_oldplasma
    use mod_sl_input
    IMPLICIT none

    ! jdemod - this entire routine is being rewritten to support dynamically allocated arrays and completely
    !          unstructured input files. Any allocated array is initialized to 0 at allocation by the allocate_arrays  
    !          routines. Only initializations to non-zero values will need to be considered.



    INTEGER i1,i2,ir


    COMMON /OPTTEMP/ osm_matcht,forcet1
    INTEGER          osm_matcht,forcet1

    !COMMON /PEAKSHAPECOM/ ipeakshape
    !INTEGER               ipeakshape(2)

    !COMMON /OSMSHAPE/ osmsplim   ,osmpksym
    !REAL              osmsplim(2),osmpksym(2)

    !      COMMON /OSMOPTS/ osm_fixopt
    !      INTEGER          osm_fixopt(2,MAXNRS)


    COMMON /TIMESTATS/ timeint
    REAL               timeint(100)

    !      COMMON /OSMPMKCOM/ save_osmpmk,scaleplateau
    !      REAL*8             save_osmpmk(0:MAXNKS,MAXNRS)
    !      LOGICAL            scaleplateau(MAXNRS)

    ntaglist = 0
    DO i1 = 1, MAXTAG
       taglist(i1) = '   '
    ENDDO

    WRITE(machine2,'(64X)')

    CALL LoadEIRENEAtomicData

    ! Initialization to zero is not required for dynamic arrays since
    ! the allocation routines take care of setting them all to zero

    !CALL ISet(nvertp,MAXNKS*MAXNRS,0)
    !CALL RZero(flxhw7,MAXSEG)
    timeint = 0.0
    !CALL RZero(timeint,100)
    !CALL IZero(virtag,MAXNKS*MAXNRS)

    sldebug = 0
    connected = .FALSE.

    rvertp = 0.0
    zvertp = 0.0

    idring = TARTOTAR
    irorg2 = 0
    ikmidplane = 0
    psindat = 0
    psidat = 0.0
    ! jdemod - change default PIN code from NIMBUS=0 to EIRENE07=5
    ! Need to verify other unstructured defaults are set appropriately
    ! Note: pincode should also define the default script to run EIRENE
    ! pincode = 0

    pincode = 5

    disindex = -1
    !CALL ISet(disindex,MAXNRS,-1)
    !CALL RZero(dislev,MAXNRS)

    indasd = 0
    !CALL RZero(pinasd,MAXASCDAT*MAXASD2*MAXASS*2)
    !...  Set additional cell "plasma" values to EIRENE vacuum defaults:
    DO i1 = 1, MAXASS
       DO i2 = 1, MAXASCDAT
          pinasd(i2,1,i1,1) = 1.0E+01
          pinasd(i2,2,i1,1) = 1.0E-02
          pinasd(i2,3,i1,1) = 0.0
          pinasd(i2,4,i1,1) = 0.0
          pinasd(i2,5,i1,1) = 0.0
       ENDDO
    ENDDO

    asc_ncell = 0
    !CALL IZero(asc_cell  ,MAXASC)
    !CALL IZero(asc_region,MAXASC)
    !CALL IZero(asc_link  ,MAXASC*4)
    !CALL IZero(asc_grid  ,MAXASC*2)
    !CALL IZero(asc_nvp   ,MAXASC)
    !CALL RZero(asc_rvp   ,MAXASC*8)
    !CALL RZero(asc_zvp   ,MAXASC*8)
    !CALL RZero(ascdata   ,MAXASCDAT*5)

    idring = TARTOTAR
    !CALL ISet(idring,MAXNRS,TARTOTAR)

    !CALL IZero(supflx,MAXNRS*2)

    mulrec = 1.0
    mulrecl= 1.0
    mulion = 1.0
    mulqer = 1.0
    mulqei = 1.0
    !CALL RSet(mulrec,MAXNKS*MAXNRS,1.0)
    !CALL RSet(mulrecl,MAXNRS,1.0)
    !CALL RSet(mulion,MAXNKS*MAXNRS,1.0)
    !CALL RSet(mulqer,MAXNKS*MAXNRS,1.0)
    !CALL RSet(mulqei,MAXNKS*MAXNRS,1.0)

    !CALL IZero(osm_ncnt ,MAXNRS)
    !CALL IZero(osm_nerr ,MAXNRS*2)
    !CALL IZero(osm_cadj ,MAXNRS*2*5)
    !CALL IZero(osm_tcon ,MAXNRS)
    !CALL RZero(osm_temod,MAXNRS)
    !CALL RZero(osm_relfr,MAXNRS)

    !CALL IZero(ikbreak,MAXNRS)

    ! logical variables are initialized to false when allocated
    !pininit = .false.
    !DO i1 = 1, MAXNRS
    !  pininit(i1) = .FALSE.
    !ENDDO

    WRITE(50,*) 'Initialising unstructured input options: '
    WRITE(50,*) ' '


    cgridst  = 0

    !CALL IZero(eiraout,MAXASD*10) 

    eirscale = 1.0
    !CALL RSet(eirscale,100,1.0)

    setmach = 0.0
    muldensity = 1.0
    addtemp = 99.0

    eirntrans = 0
    eirniontime = 0

    !...  Remove from common block.  Replaced with EIRNTORSEG.
    eirnttra = 0

    eirnaout = 1
    eiraout(1,1) =  2
    eiraout(1,2) = -MAXASCDAT
    eiraout(1,3) =  MAXASCDAT
    eiraout(1,4) =  0
    eiraout(1,5) =  0

    eirntally  =  0
    DO i1 = 1, MAXTALLY
       DO i2 = 1, 10
          WRITE(eirtally(i1,i2),'(72X)')
       ENDDO
    ENDDO
    eirphoton  =  0
    ! jdemod - these default values REQUIRE an option in the input 
    ! file since they are not valid defaults - which breaks the concept
    ! for unstructured inputs not being required
    !eirdata    = -1
    !eirgeom    = -1
    !eirgrid    = -1
    !eiradd     = -1
    !eirneut    = -1
    ! set default run time to 30 seconds for testing
    eirtime    =  30
    eirdata    =  1 
    eirgeom    =  1
    eirgrid    =  0
    ! EIRADD doesn't appear to do anything
    eiradd     = -1
    ! EIRNEUT=1 - seamless wall option (?)
    eirneut    =  1     

    !eirtrim    =  0
    ! Turn on trim database for default EIRENE run
    eirtrim    =  1

    !eirniter   =  1
    eirniter   =  0

    ! These eirmat default values are invalid - min for input is 1
    ! 2=carbon - use as default. 1-Mo 2-C 3-W 4-Be
    !eirmat1    =  0
    !eirmat2    =  0
    eirmat1    =  2
    eirmat2    =  2
    eirdebug   =  0
    eirtemp1   =  0.0
    eirtemp2   =  0.0
    eirnpgdat  =  0
    eirnasdat  =  0
    eirnspdat  =  0
    eirspdatmode = 0
    eirnpuff   =  0
    eirbgk     =  0
    eiropacity =  0
    eircxd2    =  0
    eiracx     =  1
    eirpmode   =  0
    !... Should really be 1?
    eird2ion = 1
    eird2dis = 1
    eirph2 = 0
    eiralloc   = 1.0
    eirnstrata = 0
    eirnsection = 1
    CALL IZero(eirnlimi,5000)
    CALL RZero(eirstrata,MAXSTRATA*10)
    eirntorseg = 0
    eirermin   =  1.0
    eirrinteg  =  0.0
    eireinteg  =  0.0
    eirtrc1    =  0
    eirtrc2    =  0
    eird2fu    =  0
    eirzaa     = -1.0
    eirtorfrac =  1.0
    eirsrcmul  =  1.0
    eirfuji    =  0
    eirnsdtor  = 1
    !CALL RZero(eirsdtor,MAXTOR)


    !...  Time-dependent mode parameters:
    eirdtimv = 0.0

    vacnseg   = 0
    !CALL IZero(vacregion,MAXVACREGION)

    !CALL RZero(eirpgdat,MAXNAS*MAXASD)
    !CALL RZero(eirasdat,MAXNAS*MAXASD)
    !CALL RZero(eirspdat,MAXNAS*MAXASD)

    eirnres = 0
    !      CALL RZero(eirres,6*7*MAXPINITER)
    !      CALL RZero(eirres,20*MAXPINITER)

    cmodopt  = 0
    stagopt  = 0

    thesis = .FALSE.

    stopopt  = 0
    stopopt2 = 0
    stopopt3 = 0

    !CALL IZero(iflexopt,MAXFLEX)
    !CALL RZero(rflexopt,MAXFLEX)
    !...Presently, allowed over/underionisation in AdjustParticleSource,
    !   i.e. +-10%:
    rflexopt(1) = -99.0
    rflexopt(2) = 0.00
    !      rflexopt(2) = 0.50
    rflexopt(3) = 1.00
    rflexopt(4) = 1.00
    rflexopt(5) = 1.00

    outwall   = 1
    outplasma = 1
    outcell   = 1
    outpoly   = 1
    outgeom   = 1
    outtarget = 1
    outpin    = 1
    outmode   = 0

    out_source = 0
    out_plasma = 0
    out_geom   = 0

    tarninter = NULL
    tarinter = 0.0
    tarintermode = 0.0

    psitarg = 0.0

    tarsource  = 0
    tarshift(IKLO) = 0.0
    tarshift(IKHI) = 0.0

    haldata    = 0

    grd_minpl = 1.0E-4
    grd_refine = 0
    grd_range = 1.0
    grd_shift = 0.0

    firstshape = .TRUE.
    !ipeakshape(IKLO) = MAXNKS
    !ipeakshape(IKHI) = MAXNKS

    osmbulkn = 0.0
    osmbulkte = 0.0
    osmbulkti = 0.0
    osmbulkv = 0.0


    !osmsplim(IKLO) = 0.0
    !osmsplim(IKHI) = 0.0
    !osmpksym(IKLO) = 0.0
    !osmpksym(IKHI) = 0.0

    osm_range  = 1.0
    osm_store  = -1
    osm_watch1 = 0
    osm_watch2 = 0
    osm_probe  = 0
    osm_nflag  = 0
    osm_mode   = 0
    osm_matchs = 0
    osm_matchp = 0
    osm_preopt = 0                     
    osm_recopt = 0

    !CALL IZero(osmikshift,2*MAXNRS)

    osm_testep = 0.0

    s21_mode   = 0

    s28mode = 0.0
    s28ion = 0
    s28rec = 1
    s28cfp = 0
    s28cfpnt = 0
    s28cfpdrft = 0
    s28ionpfz = 0
    s28recpfz = 1
    s28cfppfz = 1
    s28cfppfznt = 0
    s28cfppfzdrft = 0
    s28mom = 0
    s28cfm = 0
    s28mompfz = 0
    s28cfmpfz = 0
    s28te = 0
    s28ti = 0
    s28tipfz = 0
    s28tiratio = 1.0
    s28fit = 0
    s28superdet = 1
    s28superdetpfz = 1
    s28probe = 0
    s28probepfz = 0
    s28nemode = 0
    s28nemodepfz = 0
    s28momfr = 0.3
    s22pfzion = 0
    s28b34 = 2.0 / 4.0
    s28b34det = 3.0 / 4.0
    s28momTe = 10.0
    s28sym = 1

    s28recfrac = 1.0
    s28ionfrac = 1.0
    s28momfrac = 1.0

    s28ionset = .FALSE.
    s28recset = .FALSE.
    s28momset = .FALSE.
    s28ionsetpfz = .FALSE.
    s28recsetpfz = .FALSE.
    s28momsetpfz = .FALSE.

    !...  slcom_sol28:
    new = .FALSE.

    rel_iter   =  1
    rel_count  = -1
    rel_opt    =  0
    relmode    =  0
    rel_frac   =  1.0
    rel_nstep  =  0
    rel_step   =  0
    rel_niter  =  0
    rel_pace   =  1.0
    rel_ndata  =  0
    rel_bound1 =  0.0
    rel_bound2 =  1.0
    relexpdat  =  0.0
    osm_symopt =  0
    rel_tol    =  1.0
    osmmock    =  0

    prb_align = 0
    prb_shift = 0.0

    DO i1 = 1, MAXNRS

       rel_deltati   (i1) = 0.0
       rel_hte       (i1) = 0.0
       rel_hti       (i1) = 0.0
       rel_hne       (i1) = 0.0
       rel_symfr     (i1) = 1.0
       rel_prbfr(IKLO,i1) = 1.0
       rel_prbfr(IKHI,i1) = 1.0

       ! osm_fixopt is initialized to zero when allocated
       !osm_fixopt(IKLO,i1) = 0
       !osm_fixopt(IKHI,i1) = 0

       DO i2 = 1, 3
          rel_hproe(i2,i1) = HI
          rel_hproi(i2,i1) = HI
       ENDDO

       osm_model(IKLO,i1) = 0
       osm_model(IKHI,i1) = 0

       osm_dp1(IKLO,i1) = 0.15
       osm_dp1(IKHI,i1) = 0.15

       osm_dp2(IKLO,i1) = 0.0
       osm_dp2(IKHI,i1) = 0.0

       osm_dp3(IKLO,i1) = 1.0
       osm_dp3(IKHI,i1) = 1.0

       osm_dp5(IKLO,i1) = 0.0
       osm_dp5(IKHI,i1) = 0.0

       osm_code(IKLO,i1) = 1
       osm_code(IKHI,i1) = 1

       DO i2 = 1, MAXNKS
          osm_dp4(i2,i1) = 1.0
          osm_dp6(i2,i1) = 1.0
       ENDDO

       rel_qemul (IKLO,i1) = 1.0
       rel_qemul (IKHI,i1) = 1.0
       !...temp1
       osm_peimul(IKLO,i1) = -0.1
       osm_peimul(IKHI,i1) = -0.1

       ! initialized to zero at allocation
       !cfs_mage(i1) = 0.0
       !cfs_magi(i1) = 0.0
       !cfs_sume(i1) = 0.0
       !cfs_sumi(i1) = 0.0
    ENDDO

    nlpdati2  =  0
    nlpdato2  =  0

    adp_opt    = 0
    adp_region = 3
    adp_upper  = HI
    adp_lower  = 1.0

    ! jdemod - the con_ variables appear unused anywhere
    !con_nerr = 0
    DO i1 = 1, MAXNRS
       cmachno(i1,1) = 1.0
       cmachno(i1,2) = 1.0

       !con_errc (i1) = 0
       !con_errir(i1) = 0

       idds(i1,1) = MAXNDS
       idds(i1,2) = MAXNDS

       ikerror(IKLO,i1) = MAXNKS
       ikerror(IKHI,i1) = MAXNKS
    ENDDO

    korpg = maxnks * maxnrs
    !DO i1 = 1, MAXNRS
    !  DO i2 = 1, MAXNKS
    !    korpg(i2,i1) = MAXNKS*MAXNRS
    !  ENDDO
    !ENDDO

    !CALL RZero(lpdati2,MAXINS*9)
    !CALL RZero(lpdato2,MAXINS*9)

    DO i1 = 1, MAXINS
       lpdati2(i1,8) = 1.0
       lpdati2(i1,9) = 1.0
       lpdato2(i1,8) = 1.0
       lpdato2(i1,9) = 1.0

       DO i2 = IKLO, IKHI
          rel_mfsp(i2,i1) = 1.0
       ENDDO
    ENDDO

    rel_data = 1.0
    !CALL RSet (rel_data ,MAXINS*15,1.0)
    !CALL IZero(osm_sympt,MAXNRS)
    !CALL IZero(rel_viter,(1+MAXSTEP))

    !     jdemod - These are initialized to zero when allocated
    !     CALL RZero(oldktebs,MAXNKS*MAXNRS)
    !      CALL RZero(oldktibs,MAXNKS*MAXNRS)
    !      CALL RZero(oldknbs ,MAXNKS*MAXNRS)
    !      CALL RZero(oldkvhs ,MAXNKS*MAXNRS)

    s21_ndatai = 0
    s21_ndatao = 0

    !CALL RZero(s21_datai,MAXNRS*9)
    !CALL RZero(s21_datao,MAXNRS*9)

    !CALL RZero(bgplasopt,2*MAXNRS*9)

    !      switch(SWPOW2) = 1.0
    !      switch(SWPOW3) = 4.0
    !      switch(SWION2) = 1.0

    !CALL RZero(prp_te,MAXNRS*NUMPRB)
    !CALL RZero(prp_ti,MAXNRS*NUMPRB)
    !CALL RZero(prp_ne,MAXNRS*NUMPRB)

    !CALL RZero(dip_te,MAXPRB*NUMPRB)
    !CALL RZero(dip_ti,MAXPRB*NUMPRB)
    !CALL RZero(dip_ne,MAXPRB*NUMPRB)
    !CALL RZero(dip_s ,MAXPRB*NUMPRB)
    !CALL RZero(dip_v ,MAXPRB*NUMPRB)

    !CALL IZero(prb_num,NUMPRB)

    !CALL RZero(prb_te ,MAXPRB*NUMPRB)
    !CALL RZero(prb_ti ,MAXPRB*NUMPRB)
    !CALL RZero(prb_ne ,MAXPRB*NUMPRB)
    !CALL RZero(prb_rho,MAXPRB*NUMPRB)
    !CALL RZero(prb_r  ,MAXPRB*NUMPRB)
    !CALL RZero(prb_z  ,MAXPRB*NUMPRB)

    !CALL RZero(pinior ,MAXNKS*MAXNRS)
    !CALL RZero(pinior2,MAXNKS*MAXNRS)
    !CALL RZero(pinmpr ,MAXNKS*MAXNRS)
    !CALL RZero(pinmpr2,MAXNKS*MAXNRS)
    !CALL RZero(pinqir ,MAXNKS*MAXNRS)
    !CALL RZero(pinqir2,MAXNKS*MAXNRS)
    !CALL RZero(pinqer ,MAXNKS*MAXNRS)
    !CALL RZero(pinqer2,MAXNKS*MAXNRS)

    !CALL RZero(pinion  ,MAXNKS*MAXNRS)
    !CALL RZero(pinrec  ,MAXNKS*MAXNRS)
    !CALL RZero(pinqi   ,MAXNKS*MAXNRS)
    !CALL RZero(pinqe   ,MAXNKS*MAXNRS)
    !CALL RZero(pinmol  ,MAXNKS*MAXNRS)
    !CALL RZero(pinena  ,MAXNKS*MAXNRS)
    !CALL RZero(pinenz  ,MAXNKS*MAXNRS)
    !CALL RZero(pinenm  ,MAXNKS*MAXNRS)
    !      CALL RZero(pinvdist,...)
    !CALL RZero(pinmp   ,MAXNKS*MAXNRS)
    !CALL RZero(pinatom ,MAXNKS*MAXNRS)
    !CALL RZero(pinalpha,MAXNKS*MAXNRS)
    !CALL RZero(pinline(1,1,1,H_BALPHA),MAXNKS*MAXNRS*6)
    !CALL RZero(pinline(1,1,1,H_BGAMMA),MAXNKS*MAXNRS*6)
    !CALL RZero(pinionz ,MAXNKS*MAXNRS)
    !CALL RZero(pinz0   ,MAXNKS*MAXNRS)
    !CALL RZero(pindata ,MAXNKS*MAXNRS*MAXDATA)
    !CALL RZero(pindata2,MAXNKS*MAXNRS*MAXDATA)

    !CALL RZero(pinbgk  ,MAXNKS*MAXNRS*MAXBGK*MAXTOR)
    !CALL RZero(pinbgk2 ,MAXNKS*MAXNRS*MAXBGK)

    !CALL RZero(osmpei     ,MAXNKS*MAXNRS)
    !CALL RZero(osmcde     ,MAXNKS*MAXNRS)
    !CALL RZero(osmcdi     ,MAXNKS*MAXNRS)
    !CALL RZero(osmcve     ,MAXNKS*MAXNRS)
    !CALL RZero(osmcvi     ,MAXNKS*MAXNRS)
    !CALL RZero(osmpmk (0,1),(MAXNKS+1)*MAXNRS)
    !CALL DZero(osmpmk2(0,1),(MAXNKS+1)*MAXNRS)
    !CALL DZero(save_osmpmk(0,1),(MAXNKS+1)*MAXNRS)

    ! dynamic logical arrays are already initialized to false
    !     DO ir = 1, MAXNRS
    !        scaleplateau(ir) = .FALSE.
    !      ENDDO

    !      CALL RZero(pinion2  ,MAXNKS*MAXNRS)
    !      CALL RZero(pinrec2  ,MAXNKS*MAXNRS)
    !      CALL RZero(pinqi2   ,MAXNKS*MAXNRS)
    !      CALL RZero(pinqe2   ,MAXNKS*MAXNRS)
    !      CALL RZero(pinmol2  ,MAXNKS*MAXNRS)
    !      CALL RZero(pinena2  ,MAXNKS*MAXNRS)
    !      CALL RZero(pinenz2  ,MAXNKS*MAXNRS)
    !      CALL RZero(pinenm2  ,MAXNKS*MAXNRS)
    !      CALL RZero(pinvdist2,...)
    !      CALL RZero(pinmp2   ,MAXNKS*MAXNRS)
    !      CALL RZero(pinatom2 ,MAXNKS*MAXNRS)
    !      CALL RZero(pinalpha2,MAXNKS*MAXNRS)
    !      CALL RZero(pinionz2 ,MAXNKS*MAXNRS)
    !      CALL RZero(pinz02   ,MAXNKS*MAXNRS)

    !CALL RZero(osmpei2,MAXNKS*MAXNRS)

    !CALL RZero(osmcfp ,MAXNKS*MAXNRS)
    ! This is also a dynamic array - initialized at creation
    !osmcfpflx = 0.0
    !CALL RZero(osmcfi ,MAXNKS*MAXNRS)
    !CALL RZero(osmcfe ,MAXNKS*MAXNRS)

    !CALL RZero(osmmp  ,MAXNKS*MAXNRS)
    !CALL RZero(osmmp2 ,MAXNKS*MAXNRS)
    !CALL RZero(osmqe  ,MAXNKS*MAXNRS)
    !CALL RZero(osmqe2 ,MAXNKS*MAXNRS)

    gtarg = 1.0
    elecptarg = 1.0
    ionptarg = 1.0

    !DO i1 = 1, MXSPTS
    !  DO i2 = 1, 3
    !    gtarg    (i1,i2) = 1.0
    !    elecptarg(i1,i2) = 1.0
    !    ionptarg (i1,i2) = 1.0
    !  ENDDO
    !ENDDO

    !inputflag = 0
    !DO i1 = 1, 100
    !  inputflag(i1) = 0
    !ENDDO

    !...  Tags:
    tagmulrec = .FALSE.
    tagpinatom = .FALSE.

    !
    !     Call routine to initialize the unstructured input data 
    !     This is a separate routine since unlike normal variable 
    !     initialization - the values set in this routine can have
    !     a significant impact on the simulation and thus the default
    !     values need to well documented.
    !
    call InitializeUnstructuredInput

    CALL divInitializeOptions
    CALL osm_InitializeOptions
    
    return
  END SUBROUTINE InitializeVariables




  ! ====================================================================== 

  ! subroutine: InitializeUnstructuredInput 

  subroutine InitializeUnstructuredInput 
    use subgrid_options 
    use ribbon_grid_options 
    use mod_sol22_input 
    use ero_interface 
    use allocatable_input_data 
    ! slmod begin 
    use mod_grid 
    ! slmod end 
    use mod_params 
    use mod_slcom 
    use mod_cgeom 
    use mod_cadas 
    use mod_comtor 
    use mod_cedge2d 
    use mod_solparams 
    use mod_solswitch 
    use mod_solcommon 
    use mod_reiser_com 
    use mod_line_profile 
    use mod_driftvel 
    use mod_diagvel 
    use mod_fperiph_com 
    use mod_dperpz 
    use mod_lambda 
    use mod_sol23_input
    use mod_sol29_input 
    use mod_promptdep 
    use mod_rundiv_local
    use mod_pindata
    use mod_adpak_com
    use mod_dynam4
    use mod_dynam5
    implicit none 

    !     INCLUDE 'params' 

    !     INCLUDE 'slcom' 
    !     INCLUDE 'cgeom' 
    !     include 'cadas' 
    !     INCLUDE 'comtor' 
    !      INCLUDE 'pindata' 
    !     INCLUDE 'cedge2d' 

    !     INCLUDE 'solparams' 
    !     INCLUDE 'solswitch' 
    !     INCLUDE 'solcommon' 

    !     include 'reiser_com' 
    !     include 'line_profile' 

    !     include 'driftvel' 
    !     include 'fperiph_com' 
    !     include 'dperpz' 
    !      include 'slcom_sol28' 

    integer :: in


    !     Intializing Unstructured input data to default values. 


    call sol22_initialize_unstructured_input 
    call sol23_initialize_unstructured_input
    call sol29_initialize_unstructured_input 




    !call initialize_tag_series_A
    !call initialize_tag_series_B
    !call initialize_tag_series_C
    !call initialize_tag_series_D
    !call initialize_tag_series_F
    !call initialize_tag_series_G
    !call initialize_tag_series_H
    !call initialize_tag_series_I
    !call initialize_tag_series_K
    !call initialize_tag_series_N
    !call initialize_tag_series_P
    !call initialize_tag_series_Q
    !call initialize_tag_series_R
    !call initialize_tag_series_S
    !call initialize_tag_series_T
    !call initialize_tag_series_W
    !call initialize_tag_series_0


    ! slmod begin 
    !     TAG 077 

    !     Data for additional neutral wall surfaces 

    n_grid_wall = 0  ! related, so sticking here for now... -SL, 15/11/2013 

    eirnasdat = 1 
    eirasdat  = 0.0 
    eirasdat(1,1) = 998.0  ! This triggers the automated wall clipping for EIRENE 
    ! slmod end 

    !     End of intialization 



    !      subroutine initialize_tag_series_A
    !      implicit none




    !
    !  Initialization for TAG series A
    !
    !
    ! TAG:  +A01    Initialization
    !
    ! Sample Input
    !
    !  '+A01 Ref DIV Title' 'Case name and header'
    !
    ! Input call documentation
    !
    !   'TITLE FOR RUN'
    !
    TITLE = 'SPECIFY TAG A01: Case name and header are missing'
    !
    ! TAG:  +A02    Initialization
    !
    ! Sample Input
    !
    !  '+A02 DIV Desc:' 'Text Describing the case'
    !
    ! Input call documentation
    !
    !   'DESCRIPTION OF RUN'
    !
    desc = 'SPECIFY TAG A02: Text describing the case is missing'
    !
    ! TAG:  +A03    Initialization
    !
    ! Sample Input
    !
    !  '+A03 Equil File Name' 'sonnet_oskn_c08'
    !
    ! Input call documentation: NOTE: At the present time this input is documentation only since the grid file is specified on the rundiv command line
    !
    !  'equilibrium file name'
    !
    equil = 'SPECIFY TAG A03: Name of grid file used for case (DOCUMENTATION ONLY)'
    !
    ! TAG:  +A04    Initialization
    !
    ! Sample Input
    !
    !  '+A04    Print option  (0 reduced, 1 full)               '      0
    !
    ! Input call documentation
    !
    !  'PRINT OPTION'
    !
    CPRINT = 0

    !----------------------------------------------------------------------- 
    !     TAG A05 
    !     ne_opt - selects how electron density is defined 
    !     0 -> ne = nb    1 -> ne = nb + sigma nz (from fluid code) 

    ne_opt   = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG A06 

    !     Option to write a JET TRAN file for POST-PROCESSOR use 
    !     Written to Unit 41 - 0 = off - 1 = 0n .. default OFF 

    write_tran = 0 


    ! ----------------------------------------------------------------------- 

    !     TAG A07: 


    !     Option to write a netcdf version of the raw data file 
    !     0 = off   1 = 0n .. default ON 


    netcdf_opt = 1 

    !  end subroutine initialize_tag_series_A

    !    subroutine initialize_tag_series_B
    !    implicit none

    !
    !  Initialization for TAG series B
    !
    !
    ! TAG:  +B01    Initialization
    !
    ! Sample Input
    !
    !  '+B01    Special plasma parameter           Rspec        '     1
    !
    ! Input call documentation
    !
    !  'SPEC PLASMA RSPEC'
    !
    IRSPEC = 1
    !
    ! TAG:  +B02    Initialization
    !
    ! Sample Input
    !
    !  '+B02    For average density "near" target  xnear (m)    '    0.6
    !
    ! Input call documentation
    !
    !  'NEAR PARAM  XNEAR '
    !
    CXNEAR = 0.6
    !
    ! TAG:  +B03    Initialization
    !
    ! Sample Input
    !
    !  '+B03    For average density "near" target  ynear (m) +/-'    0.6
    !
    ! Input call documentation
    !
    !  'NEAR PARAM  YNEAR '
    !
    CYNEAR = 0.6
    !
    ! TAG:  +B04    Initialization
    !
    ! Sample Input
    !
    !  '+B04    Debug atoms (0 off, >0 print every nth timestep)'      0
    !
    ! Input call documentation
    !
    !  '*** DEBUG NEUT ***'
    !
    ISTEPN  = 0
    !
    ! TAG:  +B05    Initialization
    !
    ! Sample Input
    !
    !  '+B05    Debug ions  (0 off, >0 print every nth timestep)'      0
    !
    ! Input call documentation
    !
    !  '*** DEBUG DIV  ***'
    !
    ISTEPL  = 0
    !
    ! TAG:  +B06    Initialization
    !
    ! Sample Input
    !
    !  '+B06    Debug ion velocity  (0 off, >1 on)              '      0
    !
    ! Input call documentation
    !
    !  '**DEBUG VELOCITY**'
    !
    ISTEPV  = 0
    !
    ! TAG:  +B07    Initialization
    !
    ! Sample Input
    !
    !  '+B07 TN442 Z-value defining divertor region  (m)         '   1.7
    !
    ! Input call documentation
    !
    !   'Z DIVERTOR LIMIT'
    !
    CZD    = 1.7
    !
    ! TAG:  +B08    Initialization
    !
    ! Sample Input
    !
    !  '+B08 TN    Ring number for detailed background data      '    0
    !
    ! Input call documentation
    !
    !  'HI-RES RING     '
    !
    CIRHR = 0

    ! end subroutine initialize_tag_series_B



    !    subroutine initialize_tag_series_C
    !      implicit none
    !
    !  Initialization for TAG series C
    !
    !
    ! TAG:  +C01    Initialization
    !
    ! Comments from input
    !
    ! Example
    !'TN982     Number of S-values  :-'  4
    !   5.0
    !   10.0
    !   15.0
    !   20.0
    !
    ! Sample Input
    !
    !  '+C01 ' 'Set of S-distances for ion leakage diagnostic(m)'
    !  'TN982     Number of S-values  :-'  0
    !
    ! Input call documentation
    !
    !  'LEAK- S'
    !
    cleaksn = 0

    !
    ! TAG:  +C02    Initialization
    !
    ! Sample Input
    !
    !  '+C02 TN1303 Dperp Extractor - methods used - 0 1 2       '     2
    !
    ! Input call documentation
    !
    !  'DPERP EXT METHOD'
    !
    dpmethod = 2
    !
    ! TAG:  +C03    Initialization
    !
    ! Sample Input
    !
    !  '+C03 TN1310 Dperp Ext. 0=only to Xpoint 1=Full Field Line'     1
    !
    ! Input call documentation
    !
    !  'DPERP EXT SUM LIMIT'
    !
    dpsuml = 1
    !
    ! TAG:  +C04    Initialization
    !
    ! Sample Input
    !
    !  '+C04 TN1310 Dperp Ext. Outer Ring Losses    0=off 1=on   '     1
    !
    ! Input call documentation
    !
    !  'DP EXT OUTER RING'
    !
    dpouter = 1
    !
    ! TAG:  +C05    Initialization
    !
    ! Sample Input
    !
    !  '+C05 TN1309 Dperp Ext. Dperp  Convection    0=off 1=on   '     1
    !
    ! Input call documentation
    !
    !  'DP EXT CONVECT LOSS'
    !
    dpconv = 1
    !
    ! TAG:  +C06    Initialization
    !
    ! Sample Input
    !
    !  '+C06 TN1311 Dperp Ext. 1/2 Cell Flux Corr.  0=off 1=on   '     0
    !
    ! Input call documentation
    !
    !  'CELL CENTRE FLUX'
    !
    dpfluxopt = 0
    !
    ! TAG:  +C07    Initialization
    !
    ! Sample Input
    !
    !  '+C07 TN1311 Dperp Ext. Calc. Average Dperp. 0=off N=on   '     1
    !
    ! Input call documentation
    !
    !  'AVERAGE DPERP OPT'
    !
    dpavopt = 1
    !
    ! TAG:  +C08    Initialization
    !
    ! Sample Input
    !
    !  '+C08 TN1311 Dperp Ext. Major Radius Corr.   0=off 1=on   '     0
    !
    ! Input call documentation
    !
    !  'MAJOR RADIUS CORR'
    !
    dprcopt = 0
    !
    ! TAG:  +C09    Initialization
    !
    ! Sample Input
    !
    !  '+C09 TN1314 Dperp Ext. Gradient Smoothing   0=off N=on   '     0
    !
    ! Input call documentation
    !
    !  'GRADIENT SMOOTHING'
    !
    dpsmooth = 0
    !
    ! TAG:  +C10    Initialization
    !
    ! Sample Input
    !
    !  '+C10 TN     Dperp Ext. Gradient Calc Meth     -1,0,1     '     0
    !
    ! Input call documentation
    !
    !     'GRADIENT CALC METH'
    !
    dpnav = 0
    !
    ! TAG:  +C11    Initialization
    !
    ! Sample Input
    !
    !  '+C11 TN     Dperp Ext. Cross-field Area  0-centre 1-bound'     0
    !
    ! Input call documentation
    !
    !  'CELL BOUND AREA'
    !
    dparea = 0
    !
    ! TAG:  +C12    Initialization
    !
    ! Sample Input
    !
    !  '+C12 TN     Dperp Ext. Power Loss Terms     0-off 1-on   '     1
    !
    ! Input call documentation
    !
    !  'POWER LOSS TERMS'
    !
    dpploss = 1
    !
    ! TAG:  +C13    Initialization
    !
    ! Sample Input
    !
    !  '+C13 TN     Dperp Ext. Non-ortho Correction 0-off 1-on   '     1
    !
    ! Input call documentation
    !
    !  'DP NON-ORTH GRADIENT'
    !
    dporth = 1
    !
    ! TAG:  +C14    Initialization
    !
    ! Sample Input
    !
    !  '+C14 TN     Dperp Ext. Pei Correction Factor             '    1.0
    !
    ! Input call documentation
    !
    !  'PEI Correction'
    !
    dppei = 1.0
    !
    ! TAG:  +C15    Initialization
    !
    ! Sample Input
    !
    !  '+C15 TN1373 Dperp Ext. Source Recycle Correction Factor  '    1.0
    !
    ! Input call documentation
    !
    !  'EXT Recycle Frac'
    !
    dprec = 1.0
    !
    ! TAG:  +C16    Initialization
    !
    ! Sample Input
    !
    !  '+C16 TN1445 Dperp Ext. Dperp/Xperp Fixed Ratio  0.0=off  '    0.0
    !
    ! Input call documentation
    !
    !  'Dp/Xp Ratio'
    !
    dpxpratio = 0.0
    !
    ! TAG:  +C17    Initialization
    !
    ! Comments from input
    !
    !
    !     Reciprocating/Fast Scanning Probe - R location.
    !
    !     .... and horizontal probe location
    !
    !
    ! Sample Input
    !
    !  '+C17 Vertical   Reciprocating Probe - Intersection #  '     1
    !
    ! Input call documentation
    !
    !   'Crossing number'
    !
    rlocnum = 1
    !
    ! TAG:  +C18    Initialization
    !
    ! Sample Input
    !
    !  '+C18 Vertical   Reciprocating Probe - R-Value         '    1.94
    !
    ! Input call documentation
    !
    !  'Rec.probe loc'
    !
    crploc = 1.94
    !
    ! TAG:  +C19    Initialization
    !
    ! Sample Input
    !
    !  '+C19   Horizontal Reciprocating Probe - Intersection #  '     1
    !
    ! Input call documentation
    !
    !   'Crossing number'
    !
    zlocnum = 2
    !
    ! TAG:  +C20    Initialization
    !
    ! Sample Input
    !
    !  '+C20   Horizontal Reciprocating Probe - Z-Value         '    0.0
    !
    ! Input call documentation
    !
    !  'Rec.probe loc'
    !
    ! Default values are set for the location of the DIIID midplane reciprocating probe
    !
    czploc = -0.04199


    ! ----------------------------------------------------------------------- 

    !     TAG C21: 

    !     Initialization of PINQE multiplier for the Dperp extractor - this 
    !     value is usually 1.0 - however, by allowing for this multiplier 
    !     it is possible to adjust for impurity radiated power when this 
    !     information is not available. 

    dp_pinqe_mult = 1.0 

    ! ----------------------------------------------------------------------- 

    !     TAG C22 

    !     line_profile_opt is set to 0 (off) by default. There are a number 
    !     of input values required when this option is active and these 
    !     input values MUST immediately follow the unstructured input option. 

    !     The first line contains a set of ADAS selectors in the standard 
    !     format. 

    !     'Text'   'ADASID'  ADASYR 'ADAS EXTENSION'  ISELE ISELR ISELX ISELD 
    !     - only ISELE is used at the present time. 

    !     The second input line defines the LOS for the calculation and 
    !     the instrument and bin characteristics. 

    !     'Text'  ROBS  ZOBS  THETA  DTHETA INSTRUMENT_WIDTH BIN_WIDTH 

    !     No default values are assigned to any sub-options when this 
    !     option is off - all options MUST be specified when it is on. 

    line_profile_opt = 0 

    !  end subroutine initialize_tag_series_C



    !    subroutine initialize_tag_series_D
    !      implicit none



    !
    !  Initialization for TAG series D
    !
    !
    ! TAG:  +D01    Initialization
    !
    ! Sample Input
    !
    !  '+D01    Source Data Option 0-Nocorona 1-Adas            '     1
    !
    ! Input call documentation
    !
    !   'Rad/ioniz data source'
    !
    cdatopt = 1
    !
    ! TAG:  +D02    Initialization
    !
    ! Sample Input
    !
    !  '+D02    USERID for ADAS H data (*=use central database) '   '*'
    !
    ! Input call documentation
    !
    !  'ADAS H userid'
    !
    useridh = '*'
    !
    ! TAG:  +D03    Initialization
    !
    ! Sample Input
    !
    !  '+D03    Year for ADAS H data                            '    96
    !
    ! Input call documentation
    !
    !  'ADAS H year          '
    !
    iyearh = 96
    !
    ! TAG:  +D04    Initialization
    !
    ! Sample Input
    !
    !  '+D04    USERID for ADAS Z data (*=use central database) '   '*'
    !
    ! Input call documentation
    !
    !  'ADAS Z userid'
    !
    useridz = '*'
    !
    ! TAG:  +D05    Initialization
    !
    ! Sample Input
    !
    !  '+D05    Year for ADAS Z data                            '    96
    !
    ! Input call documentation
    !
    !  'ADAS Z year          '
    !
    iyearz = 96
    !
    ! TAG:  +D06    Initialization
    !
    ! Sample Input
    !
    !  '+D06 B2FRATES name:' '/home/david/divimp/adpak/C_rates.strahl'
    !
    ! Input call documentation
    !
    !  'Name of File for ADPAK/INEL atomic data'
    !
    mcfile = '/home/david/divimp/adpak/C_rates.strahl'
    !
    ! TAG:  +D07    Initialization
    !
    ! Sample Input
    !
    !  '+D07    Sputter data option 1-old 2-93                  '     6
    !
    ! Input call documentation
    !
    !  'SPUTTER SOURCE OPT '
    !
    CSPUTOPT = 6
    !
    ! TAG:  +D08    Initialization
    !
    ! Sample Input
    !
    !  '+D08    Chemical Sputter Data Option                    '    11
    !
    ! Input call documentation
    !
    !  'CHEMSPUT SOURCE OPT'
    !
    CCHEMOPT = 11
    !
    ! TAG:  +D09    Initialization
    !
    ! Sample Input
    !
    !  '+D09 TN1481 BG Ion Mom.Tran.Coll.Coef.     (kelighi)     '    5.0E-16
    !
    ! Input call documentation
    !
    !  'MTC COEFF 1 Imp-I'
    !
    kelighi = 5.0E-16
    !
    ! TAG:  +D10    Initialization
    !
    ! Sample Input
    !
    !  '+D10 TN1481 BG Neutral Mom.Tran.Coll.Coef. (kelighg)     '    5.0E-16
    !
    ! Input call documentation
    !
    !  'MTC COEFF 2 Imp-N'
    !
    kelighg = 5.0E-16
    !
    ! TAG:  +D11    Initialization
    !
    ! Sample Input
    !
    !  '+D11    Characteristic energy Ebd          Ebd   (eV)   '    0.0
    !
    ! Input call documentation
    !
    !  'CHAR ENERGY EBD'
    !
    CEBD = 0.0
    !
    ! TAG:  +D12    Initialization
    !
    ! Sample Input
    !
    !  '+D12    Neutral hydrogen density parameter Nhc   (m**-3)'    1.0E15
    !
    ! Input call documentation
    !
    !  'H DENSITY:  NHC   '
    !
    CNHC = 1.0E15
    !
    ! TAG:  +D13    Initialization
    !
    ! Sample Input
    !
    !  '+D13                                       Nho   (m**-3)'    3.0E18
    !
    ! Input call documentation
    !
    !  'H DENSITY:  NHO   '
    !
    CNHO = 3.0E18
    !
    ! TAG:  +D14    Initialization
    !
    ! Sample Input
    !
    !  '+D14                                       lamhx (m)    '    0.02
    !
    ! Input call documentation
    !
    !  'H DENSITY:  LAMHX '
    !
    CLAMHX = 0.02
    !
    ! TAG:  +D15    Initialization
    !
    ! Sample Input
    !
    !  '+D15                                       lamhy (m)    '    0.11
    !
    ! Input call documentation
    !
    !  'H DENSITY:  LAMHY '
    !
    CLAMHY = 0.11
    !
    ! TAG:  +D16    Initialization
    !
    ! Sample Input
    !
    !  '+D16    Constant for CX Recomb option 2    Vcx   (m/s)  '    2.4E4
    !
    ! Input call documentation
    !
    !  'CXREC CONSTANT VCX'
    !
    CVCX   = 2.4E4
    !
    ! TAG:  +D17    Initialization
    !
    ! Sample Input
    !
    !  '+D17    Threshold yield for self-sputtering      (eV)   '    0.1
    !
    ! Input call documentation
    !
    !  'SELF-SPU THRESHOLD'
    !
    CTRESH = 0.1
    !
    ! TAG:  +D18    Initialization
    !
    ! Sample Input
    !
    !  '+D18    Bombarding ion charge state        Zimp         '      2
    !
    ! Input call documentation
    !
    !  'BOMBION CHARGE    '
    !
    CBOMBZ = 2
    !
    ! TAG:  +D19    Initialization
    !
    ! Sample Input
    !
    !  '+D19    Bombion type 0Zb 1H 2D 3T 4He4 5C 6Zi 7O        '      5
    !
    ! Input call documentation
    !
    !  'BOMBION FLAG 0:7  '
    !
    CBOMBF = 5
    !
    ! TAG:  +D20    Initialization
    !
    ! Sample Input
    !
    !  '+D20    Ionisation rate factor for neutrals          IRF'    1.0
    !
    ! Input call documentation
    !
    !  'IONISE RATE FACTOR'
    !
    CIRF   = 1.0
    !
    ! TAG:  +D21    Initialization
    !
    ! Sample Input
    !
    !  '+D21    Sputtering Enhancement Factor                SEF'    1.0
    !
    ! Input call documentation
    !
    !  'SPUT ENHANCE FACT.'
    !
    CSEF   = 1.0
    !
    ! TAG:  +D22    Initialization
    !
    ! Comments from input
    !
    ! IF nymfs is not specified OR if the first cymfs entry doesn't have an index of 0 then 
    ! the default yield modifier is 1.0 and reflection modifier (last entry) is 0.0
    !
    !         0.0   0.0    1.0     1.0    1.0    1.0    1.0     1.0
    !
    !---- READ IN YIELD MODIFIER FUNCTIONS
    !
    ! note default values should be set so that only modifications need to be read in
    !
    ! Sample Input
    !
    !  '+D22 ' 'Set of Yield Modifiers for Primary, Secondary neutrals'
    !  '      Number of rows of (X,Mpt,Mst,Mct,Mpw,Mcw,Refl) data-'  0
    !
    ! Input call documentation
    !
    !  'SET YIELD(P,S,CT,PW,CW) VALS'
    !
    NYMFS = 0
    
    !
    ! TAG:  +D23    Initialization
    !
    ! Sample Input
    !
    !  '+D23 TN? Fixed Yield Value for Sputter Data Option 4     '     0.001
    !
    ! Input call documentation
    !
    !  'Fixed Yield'
    !
    const_yield = 0.001
    !
    ! TAG:  +D24    Initialization
    !
    ! Sample Input
    !
    !  '+D24 TN1209 Target Temperature (K) for Chem. Sputt. Opt. '   300.0
    !
    ! Input call documentation
    !
    !  'Target Temp (K) '
    !
    ctargt  = 300.0
    !
    ! TAG:  +D25    Initialization
    !
    ! Sample Input
    !
    !  '+D25       Main Wall Temperature (K) for Chem. Sputt.   '   300.0
    !
    ! Input call documentation
    !
    !  'Wall Temp (K)   '
    !
    cwallt  = 300.0
    !
    ! TAG:  +D26    Initialization
    !
    ! Sample Input
    !
    !  '+D26 TN1450 PP Wall   Temperature (K) for Chem. Sputt.   '   300.0
    !
    ! Input call documentation
    !
    !  'PP Wall Temp (K)'
    !
    cwallp  = 300.0
    !
    ! TAG:  +D27    Initialization
    !
    ! Sample Input
    !
    !  '+D27 ' 'TN1450 Wall Temperatures (K) for specific segments'
    !  '      Number of segment ranges (Index1 Index2 Temp):'      0
    !
    ! Input call documentation
    !
    !  'WALL SEGMENT TEMPERATURES'
    !
    NWLTEMP = 0

    !
    ! TAG:  +D28    Initialization
    !
    ! Sample Input
    !
    !  '+D28    Temperature Gradient Coefficient parameter  ZO  '      4
    !
    ! Input call documentation
    !
    !  'ZO PARAMETER      '
    !
    CZO    = 4
    !
    ! TAG:  +D29    Initialization
    !
    ! Sample Input
    !
    !  '+D29      CEMAXF factor for sputtering velocities       '   1.0
    !
    ! Input call documentation
    !
    !  'EMAX-FACTOR'
    !
    CEMAXF = 1.0
    !
    ! TAG:  +D30    Initialization
    !
    ! Sample Input
    !
    !  '+D30 TN521 Impact Energy for wall launch Vel. dist (eV)  '  100.0
    !
    ! Input call documentation
    !
    !  'EIMP- WALL LAUNCH'
    !
    CEIMP = 100.0
    !
    ! TAG:  +D31    Initialization
    !
    ! Sample Input
    !
    !  '+D31 TN83? Maximum Number of sputtered generations       '   75
    !
    ! Input call documentation
    !
    !   'MAX. GENERATIONS'
    !
    CMAXGENS = 75
    !
    ! TAG:  +D32    Initialization
    !
    ! Sample Input
    !
    !  '+D32 TN1007 ABSFAC or Power - Specified - Use if > 0.0   '   0.0
    !
    ! Input call documentation
    !
    !   'ABSFAC modifier'
    !
    nabsfac = 0.0
    !
    ! TAG:  +D33    Initialization
    !
    ! Sample Input
    !
    !  '+D33  TN1200 Stgrad - Distance where Tgrad forces -> 0    '   1.0
    !
    ! Input call documentation
    !
    !   'TGRAD FORCES -> 0'
    !
    Cstgrad = 1.0
    !
    ! TAG:  +D34    Initialization
    !
    ! Sample Input
    !
    !  '+D34    H Recombination Calculation Option 4-Adas 0+-oth'     4
    !
    ! Input call documentation
    !
    !  'Recomb. Calc Opt'
    !
    crecopt = 4
    !
    ! TAG:  +D35    Initialization
    !
    ! Sample Input
    !
    !  '+D35    Recombination Limit Cut-OFF Temperature (eV)    '    0.0
    !
    ! Input call documentation
    !
    !  'Recomb. cutoff T'
    !
    treccut = 0.0
    !
    ! TAG:  +D36    Initialization
    !
    ! Comments from input
    !
    !
    !     Factor for temperature gradient force modifier option
    !
    !
    ! Sample Input
    !
    !  '+D36 TN1429 T-GRAD Force Modification - Factor Applied   '    1.0
    !
    ! Input call documentation
    !
    !  'FGRAD MOD FACT'
    !
    fgradfact = 1.0


    ! ----------------------------------------------------------------------- 

    !     TAG D37 and D38 
    !     ADAS IONIZATION AND RECOMBINATION RATE MODIFIERS 

    !     These values can be used to modify the rates read in from ADAS 
    !     The default values should always be set to 1.0 and these 
    !     should be used with care if used at all. 

    !     D37 - ADAS Ionization rate multiplier 

    adas_iz_rate_mult = 1.0 

    !     D38 - ADAS Recombination rate multiplier 

    adas_rec_rate_mult = 1.0 


    ! ----------------------------------------------------------------------- 

    !     TAG D39 : Alternate Sputter data specifier - used to select one of 
    !               several custom sputter datasets - usually based 
    !               on different impact angles 

    !     Set to normal incidence data as default 

    extra_sputter_angle = 0.0 

    ! ----------------------------------------------------------------------- 

    !     TAG D40 : Flux fraction for alternate bombarding ion sputter 
    !               calculations. Goes with CBOMBF and CBOMBZ to allow for 
    !               specification of trace impurity sputtering in cases 
    !               where hydrogenic sputtering is expected to be negligible. 

    !               Ideally this should be expanded to allow for 
    !               hydrogen+a specifiable distribution of impurity charge states 

    !               To keep current code functionality the default value for this 
    !               specifier is 1.0 

    cbomb_frac = 1.0 

    !  end subroutine initialize_tag_series_D

    !    subroutine initialize_tag_series_F
    !    implicit none



    !
    !  Initialization for TAG series F
    !
    !
    ! TAG:  +F01    Initialization
    !
    ! Sample Input
    !
    !  '+F01    Read EDGE2D BG for reference  0=No  1=Yes       '     0
    !
    ! Input call documentation
    !
    !  'READ E2D BG FOR REF  '
    !
    cre2d  = 0
    !
    ! TAG:  +F02    Initialization
    !
    ! Sample Input
    !
    !  '+F02    Use EDGE2D Target Data Option 0=Off 1=Reg 2=Flux'     0
    !
    ! Input call documentation
    !
    !  'EDGE2D TARG COND'
    !
    e2dtargopt = 0
    !
    ! TAG:  +F03    Initialization
    !
    ! Sample Input
    !
    !  '+F03    Plasma condition for missing SOL rings     CNIWA'      1
    !
    ! Input call documentation
    !
    !  'CNIWA OPTION      '
    !
    CNIWA  = 1
    !
    ! TAG:  +F04    Initialization
    !
    ! Sample Input
    !
    !  '+F04 ' 'EDGE1D/2D Deuterium drift vel. mult. factor VMF '
    !  '            Number of VMF blocks                          '      0
    !  '            Ring Range :-' -20     -30
    !  '            J0 & J1    :-'   5       5
    !  '            VMF0,1,2   :-'   1.000   1.000   1.000
    !
    ! Input call documentation
    !
    !  'VEL. MULT. FACTOR'
    !
    !  ****************   Special ********************
    !
    !    Note: the code discards the three lines of input if cnvmf=0 so init values for those are not required unless vmf blocks are
    !          included in the input file.
    !          The VMF code was particularly poorly written to require at least one block even if it would be ignored.
    !          So this input should NOT be included unless at least three lines of VMF data are also present even if CNVMF is zero in the input file

    cnvmf = 0


    !
    ! TAG:  +F05    Initialization
    !
    ! Sample Input
    !
    !  '+F05 SONNET-Number of Fluids in B2 Background Plasma File'     7
    !
    ! Input call documentation
    !
    !  'NUM.FLUIDS IN B2 BG'
    !
    nfla = 7
    !
    ! TAG:  +F06    Initialization
    !
    ! Sample Input
    !
    !  '+F06       Read Aux. Background Input File    0=off 1=on'     0
    !
    ! Input call documentation
    !
    !  'READ AUX BG FILE'
    !
    readaux = 0


    ! ----------------------------------------------------------------------- 

    !     TAG F11 

    !     This is set to 1 to indicate that a UEDGE/fluid code background has 
    !     been loaded and that the fluid code ionization data has been loaded 
    !     into PIN arrays. 

    uedge_bg = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG F12: 

    !     fc_target_calc_option - this option affects the calculation  of the 
    !     target conditions that are extracted from the background plasma 
    !     of a fluid code solution. This affects the values that are assigned 
    !     to the e2dtarg array at the point when the plasma solution is read in. 

    !     Each of the sub-options can also be set explicitly to different 
    !     values. 

    !     The fc_target_calc_option's supported are: 

    !     Option 0: EDGE2D standard 
    !               fc_v_calc_opt  = 0 
    !               fc_te_calc_opt = 1 
    !               fc_ti_calc_opt = 2 
    !               fc_ne_calc_opt = 2 
    !     Option 1: UEDGE standard 
    !               fc_v_calc_opt  = 0 
    !               fc_te_calc_opt = 0 
    !               fc_ti_calc_opt = 0 
    !               fc_ne_calc_opt = 2 
    !     Option 2: Base B2/Eirene standard 
    !               fc_v_calc_opt  = 1 
    !               fc_te_calc_opt = 2 
    !               fc_ti_calc_opt = 2 
    !               fc_ne_calc_opt = 2 
    !     Option 3: Alternate B2/Eirene 
    !               fc_v_calc_opt  = 1 
    !               fc_te_calc_opt = 1 
    !               fc_ti_calc_opt = 1 
    !               fc_ne_calc_opt = 2 

    !     F13: 

    !     fc_ne_calc_opt: 0 : Ti is value in guard cell 
    !                   : 1 : Ti is arithmetic average of guard and first cell 
    !                   : 2 : Ti is value in real cell 

    !     F14: 

    !     fc_te_calc_opt: 0 : Te is value in guard cell 
    !                   : 1 : Te is arithmetic average of guard and first cell 
    !                   : 2 : Te is value in real cell 

    !     F15: 

    !     fc_ti_calc_opt: 0 : Ti is value in guard cell 
    !                   : 1 : Ti is arithmetic average of guard and first cell 
    !                   : 2 : Ti is value in real cell 

    !     F16: 

    !     fc_v_calc_opt  : 0 : Velocity is set to value at boundary from fc 
    !                    : 1 : Velocity is set to sound speed 


    !     Set defaults to match EDGE2D target condition interpretation option 

    fc_target_calc_option = 0 
    fc_v_calc_opt  = 0 
    fc_te_calc_opt = 1 
    fc_ti_calc_opt = 2 
    fc_ne_calc_opt = 2 

    ! ----------------------------------------------------------------------- 

    !     F17: 

    !     B2 and B2.5 write a fort.31 file with an identical format - however 
    !     the parallel velocity arrays are different in the two cases - for 
    !     B2 the velocity is defined on the cell edges while for B2.5 the 
    !     velocity array in this file is defined at the cell centers. This means 
    !     that the data need to be interpreted differently in the two cases. 
    !     This option may also apply to the fluxes that are written to this 
    !     file. 

    !     option 0 = cell boundary value in file 
    !     option 1 = cell center value in file 

    fc_v_interp_opt = 0 


    ! ----------------------------------------------------------------------- 

    !     TAG F18: 

    !     SOLPS or fluid code input format option 
    !     The fort.31 file from SOLPS is used to load the background plasma 
    !     from SOLPS. However, the format of this file occasionally changes 
    !     and someone added the toroidal velocity after the poloidal and 
    !     radial velocities. This option allows the code to adjust to the 
    !     altered format. 
    !     Option 0 is the default and corresponds to the SOLPS 4.3-> 5.1 version 
    !     of the fort.31 file produced by b2plot. 
    !     Option 1 is for the SOLPS-ITER version. If they add other content 
    !     more versions may be required in the 

    !     0 = SOLPS 4.3 
    !     1 = SOLPS 5.1/ITER - change sign of parallel velocity 
    !     2 = SOLPS 5.1/ITER - do not change the sign of the parallel velocity 

    e2dformopt = 0 



    ! ----------------------------------------------------------------------- 

    !     F19: 

    !     e2dneut_select = 1 

    !     If readaux = 3 then the code tries to read a SOLPS auxiliary input 
    !     file fort.44 from a SOLPS run which contains neutral density 
    !     data for D and the impurity species in the run. 

    !     This option specifies which block of neutral data to load into the 
    !     fluid code neutral density array. 

    e2dneut_select = 1 

    !     ----------------------------------------------------------------------- 

    !     F20: 

    !     e2dion_select = 1 

    !     The e2dnzs data stored in fort.31 can contain multiple fluid species 
    !     in addition to H+. This quantity specifies an offset into this data 
    !     so that the code can start reading the correct impurity into the e2d 
    !     fluid code arrays like e2dnzs. This is required in DIVIMP so it can 
    !     include meaningful comparisons between the DIVIMP and fluid code results. 

    e2dion_select = 1 




    !  end subroutine initialize_tag_series_F


    !      subroutine initialize_tag_series_G
    !      implicit none


    !
    !  Initialization for TAG series G
    !
    !
    ! TAG:  +G01    Initialization
    !
    ! Sample Input
    !
    !  '+G01    Grid Optgion    0-JET 1-ASDEX 2-ITER             '     3
    !
    ! Input call documentation
    !
    !  'GRID OPTION          '
    !
    cgridopt = 3
    !
    ! TAG:  +G02    Initialization
    !
    ! Sample Input
    !
    !  '+G02    Non-Orthogonal Grid option 0-off 1-JET N.O.     '     3
    !
    ! Input call documentation
    !
    !  'NON-ORTHO. GRID OPT  '
    !
    northopt = 3
    !
    ! TAG:  +G03    Initialization
    !
    ! Sample Input
    !
    !  '+G03    Parallel Distance Option 0-centers 1-edges      '     1
    !
    ! Input call documentation
    !
    !     'PARALLEL DIST (S) OPT'
    !
    pdopt = 1
    !
    ! TAG:  +G04    Initialization
    !
    ! Sample Input
    !
    !  '+G04    Cross-field Distance Option 0-centres 1-edges   '     1
    !
    ! Input call documentation
    !
    !  'Proper Cross-field dist'
    !
    cfdopt = 1
    !
    ! TAG:  +G05    Initialization
    !
    ! Sample Input
    !
    !  '+G05    RZ calculation option 0-centers 1-Actual RZ     '     2
    !
    ! Input call documentation
    !
    !     'CALCULATE ACTUAL R
    !
    rzopt = 2
    !
    ! TAG:  +G06    Initialization
    !
    ! Sample Input
    !
    !  '+G06    XY Grid option 0-off 1-on                       '     0
    !
    ! Input call documentation
    !
    !  'XY GRID OPTION       '
    !
    xygrid = 0
    !
    ! TAG:  +G07    Initialization
    !
    ! Sample Input
    !
    !  '+G07    Cell Area Calculation Option 0-approx 1-polygon '     1
    !
    ! Input call documentation
    !
    !   'Cell Volumes from PGs'
    !
    cvolopt = 1
    !
    ! TAG:  +G08    Initialization
    !
    ! Sample Input
    !
    !  '+G08 T   Ion Wall Option        0 to 2                   '     2
    !
    ! Input call documentation
    !
    !  'ION WALL OPTION      '
    !
    CIONR  = 2
    !
    ! TAG:  +G09    Initialization
    !
    ! Sample Input
    !
    !  '+G09 T   Neutral Wall Option    0 to 4                   '     4
    !
    ! Input call documentation
    !
    !  'NEUTRAL WALL OPTION  '
    !
    CNEUR  = 4
    !
    ! TAG:  +G10    Initialization
    !
    ! Sample Input
    !
    !  '+G10 T   Trap Wall Option       0 to 4                   '     4
    !
    ! Input call documentation
    !
    !  'TRAP WALL OPTION     '
    !
    CTRAP  = 4
    !
    ! TAG:  +G11    Initialization
    !
    ! Sample Input
    !
    !  '+G11 T   Vessel Wall Redefinition Option (Baffle Removal)'     0
    !
    ! Input call documentation
    !
    !  'VESSEL REDEF. OPT    '
    !
    redefopt = 0
    !
    ! TAG:  +G12    Initialization
    !
    ! Sample Input
    !
    !  '+G12 T   Target Position Option 0 to 6                   '     6
    !
    ! Input call documentation
    !
    !  'TARGET POSITION OPT  '
    !
    CTARGOPT = 6
    !
    ! TAG:  +G13    Initialization
    !
    ! Sample Input
    !
    !  '+G13 T   Pre-defined geometry option -1-off 0-719 1-307  '    -1
    !
    ! Input call documentation
    !
    !  'GEOMETRY OPTION      '
    !
    CGEOOPT = -1
    !
    ! TAG:  +G14    Initialization
    !
    ! Sample Input
    !
    !  '+G14    Central Mirror Ring Location       (IR)         '     1
    !
    ! Input call documentation
    !
    !  'CORE MIRROR RING SPEC' 
    !
    ircore = 1
    !
    ! TAG:  +G15    Initialization
    !
    ! Sample Input
    !
    !  '+G15    Rectangular grid for neutrals 0calculate 99file '      0
    !
    ! Input call documentation
    !
    !  'RECT GRID CALC.   '
    !
    CRECT  = 0
    !
    ! TAG:  +G16    Initialization
    !
    ! Sample Input
    !
    !  '+G16 ' 'TN    Set of Plate coordinates                  '
    !  '    TN    Ring #, Outer R,Z   , Inner R,Z :      '       0
    !
    ! Input call documentation
    !
    !  'PLATE COORDINATES'
    !
    NPLAT = 0

    !
    ! TAG:  +G17    Initialization
    !
    ! Sample Input
    !
    !  '+G17 ' 'Wall coordinates                                '
    !  '    TN    R,Z coordinates starting at outer plate   '   0
    !
    ! Input call documentation
    !
    !  'WALL COORDINATES'
    !
    NWALL = 0

    !
    ! TAG:  +G18    Initialization
    !
    ! Sample Input
    !
    !  '+G18 ' 'Wall coordinates - PFZ                          '
    !  '    TN    R,Z coordinates                           '  0
    !
    ! Input call documentation
    !
    !  'WALL COORDINATES'
    !
    NWALL2 = 0

    !
    ! TAG:  +G19    Initialization
    !
    ! Comments from input
    !
    !
    !     NOTE: These items are now usually on the GEOM: line at the beginning
    !           of the grid file or are calculated automatically from the grid
    !     The following set of numbers specify the characteristic values of
    !     the ASDEX U style grid file to be read in. The Number of rings
    !     and elements per ring (a ring is a set of polygons in a row or
    !     indexed by the same ring number), the cutring,(which specifies the
    !     end of  rings for the core and trapped plasma), and the cutpts 1
    !     and 2 which specify the points on the rings numbered less than the
    !     cut ring where the splits for core and trapped or private plasma
    !     occur.
    !
    !
    ! Sample Input
    !
    !  '+G19 ASDEX U - GRID CHARACTERISTICS:  Number of Rings    '    26
    !
    ! Input call documentation
    !
    !  'MAXRINGS in AUG'
    !
    maxrings = 26
    !
    ! TAG:  +G20    Initialization
    !
    ! Sample Input
    !
    !  '+G20                                 Number of Knots    '    34
    !
    ! Input call documentation
    !
    !  'MAX PTS in AUG'
    !
    maxkpts = 34
    !
    ! TAG:  +G21    Initialization
    !
    ! Sample Input
    !
    !  '+G21                                 Cut ring           '     7
    !
    ! Input call documentation
    !
    !  'CUTRING in AUG'
    !
    cutring = 7
    !
    ! TAG:  +G22    Initialization
    !
    ! Sample Input
    !
    !  '+G22                                 Cut point 1        '     1
    !
    ! Input call documentation
    !
    !  'CUTPT1 in AUG'
    !
    cutpt1 = 1
    !
    ! TAG:  +G33    Initialization
    !
    ! Sample Input
    !
    !  '+G33                                 Cut point 2        '    34
    !
    ! Input call documentation
    !
    !  'CUTPT2 in AUG'
    !
    cutpt2 = 34


    ! ----------------------------------------------------------------------- 

    !     TAG G23: 

    !     This option is used to tag a SONNET style grid as being an 
    !     FRC (Field Reversed Configuration) custom grid. At the moment 
    !     only one option is supported but further subtypes could be 
    !     added for Sonnet grids requiring special processing. 

    !     0 = stanard Sonnet grid = Default 
    !     1 = FRC version 1 - type of Sonnet grid 
    !         - used to set various FRC related options 
    !     2 = sonnet grid without boundary cells - boundary cells are added 
    !         - useful for carre grids 

    sonnet_grid_sub_type = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG G34: 

    !     This input value is used to define the machine type for 
    !     the purpose of placing data in the TRAN file for JET 
    !     post-processor use. The default value is for JET. Unfortunately, 
    !     the grid option itself is not sufficiently selective since it defines 
    !     formats of grids that may or may not be associated with specific 
    !     machines. 

    !       0 = jet 
    !       1 = diiid 
    !       2 = alcator (cmod) 
    !       3 = aug (asdex upgrade) 
    !       4 = iter 
    !       5 = ignitor 
    !       6 = fire 

    tmachine_opt = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG G35: 

    !     This input value is used to define a shotid 
    !     for the case for cataloguing in the JET catalog 
    !     system. This must include the shot number. If the 
    !     shot number is present in the case name then this 
    !     value is not required - the case name may be used as 
    !     the shotid in that case. This input takes a default 
    !     value of ' '. Meaning that the shot number is expected to 
    !     be a part of the case name. 

    divshotid = ' ' 


    ! ----------------------------------------------------------------------- 

    !     TAG G36: 

    !     s_reflect_opt 

    !     This input value is used to activate parallel ion 
    !     reflection. This will insert mirrors for parallel ion transport 
    !     at ONE parallel break in the grid along each field line. This is 
    !     used to model ion transport on double null half grids. The default 
    !     value is set to off or zero. 

    s_reflect_opt = 0 

    ! ----------------------------------------------------------------------- 

    !     G37: Used by Steve for some sort of grid option - no default 
    !          values specified 

    ! ----------------------------------------------------------------------- 

    !     The following tags are related to the subgrid option for recording 
    !     more detailed data on a finer grid 

    !     G38: Base subgrid option ON/OFF 
    !     G39: R,Z dimensions of the grid region 
    !     G40: RMIN,RMAX of the subgrid region 
    !     G41: ZMIN,ZMAX of the subgrid region 

    !     Note: Although default values are assigned - these values 
    !           should be explicitly set when this option is used. 

    !     G38: Base subgrid option - OFF 

    subgrid_opt = 0 

    !     G39: Dimensions of subgrid 

    sg_rdim = 100 
    sg_zdim = 100 

    !     G40: RMIN and RMAX values of the gridded region 

    sg_rmin = 1.0 
    sg_rmax = 1.5 

    !     G41: ZMIN and ZMAX values of the gridded region 

    sg_zmin = 0.0 
    sg_zmax = 1.4 

    ! ----------------------------------------------------------------------- 

    !     Options related to ribbon grids 
    !     G42 - grid generation option - <i4> 
    !     G43 - intersection point averaging option - opt_block_av - <r4> 
    !     G44 - maximum R separation in grid generator - max_r_sep - <r4> 
    !     G45 - maximum S/Z separation in grid generator - max_s_sep - <r4> 
    !     G46 - min number of cells on ring - min_cells - <i4> 
    !     G47 - castem output identifier - <string> 
    !     G48 - min and max S for selecting intersection subset  2 x <r4> 
    !     G49 - min and max R for intersection subset grid generation 2 x <r4> 
    !     G50 - min and max S for intersection subset grid generation 2 x <r4> 
    !     G51 - length cutoff factor for ring generation <r4> default = 0.0 
    !     G52 - Cell spacing option 
    !     G53 - cell spacing factor 
    !     G54 - Input file option - RAY or CASTEM 

    !------------------------------------------------------------------------ 

    !     G42 - grid option 
    !           0 = unstructured 
    !           1 = structured 
    !           default = unstructured 

    rg_grid_opt = 0 

    !     G43 - block averaging option (removes blobs of intersection data) 

    rg_block_av = 0 

    !     G44 - maximum R separation between rows 

    rg_max_r_sep = 0.002 

    !     G45 - maximum S/Z separation between cells along row 

    rg_max_s_sep = 0.5 

    !     G46 - minimum number of cells in a row 

    rg_min_cells = 5 

    !     G47 - Castem data set to read in 

    rg_castem_data = '100610' 

    !     G48 - min and max S for selecting intersection subset  2 x <r4> 
    !           if min=max then window option is not selected 

    rg_int_win_mins = 0.0 
    rg_int_win_maxs = 0.0 

    !     G49 - min and max R for intersection subset grid generation 2 x <r4> 
    !           These are only used for subset grid generation 

    rg_minr = 0.0 
    rg_maxr = 0.0 

    !     G50 - min and max S for intersection subset grid generation 2 x <r4> 
    !           These are only used for subset grid generation 

    rg_mins = 0.0 
    rg_maxs = 0.0 

    !     G51 - length cutoff to eliminate short rings far from the separatrix from 
    !           ring generation ... some testing of this will be necessary 
    !           to obtain an optimal grid. default value is 0.0 which 
    !           effectively turns this feature off. 

    lcutoff = 0.0 

    !     G52 - Cell spacing option - default value is exponential with the exponent 
    !           specified by G53. A value of 1.0 for the cell spacing factor gives 
    !           a linear spacing. This option works better with structured grids at 
    !           the moment. 

    cell_spacing_option = 0 

    !     G53 - Cell spacing factor for determining the distribution of cells 
    !           between fixed points on rings. 
    !           default = 1.0 which gives a linear spacing 

    cell_spacing_factor = 1.0 

    !     G54 - Intersection data input file format option 
    !           0 = CASTEM 
    !           1 = RAY 

    ribbon_input_format_opt = 1 

    !------------------------------------------------------------------------ 

    !     G55 - IK offsets 

    !     jdemod - ikoffsets to move the center of the ring for 
    !              background plasma calculation 

    !     Default for all values is zero 

    n_ik_offsets = 0 
    ik_offset_data = 0.0 

    !     G56 - SOL22 half ring length option 

    !     0 - ringlen/2 
    !     1 - ksb(midnks,ir) 

    sol22_halfringlen_opt = 1 


    !  end subroutine initialize_tag_series_G

    !    subroutine initialize_tag_series_H
    !      implicit none



    !
    !  Initialization for TAG series H
    !
    !
    ! TAG:  +H01    Initialization
    !
    ! Sample Input
    !
    !  '+H01    PIN Random number seed  (<0=1, 0 generate new)  '      0
    !
    ! Input call documentation
    !
    !  'PIN RANDOM NUM. SEED'
    !
    PINISEED = 0
    !
    ! TAG:  +H02    Initialization
    !
    ! Sample Input
    !
    !  '+H02    PIN Data Print option  (0 reduced, 1 more)      '      0
    !
    ! Input call documentation
    !
    !  'PIN DATA PRINT OPT'
    !
    PINPRINT = 0
    !
    ! TAG:  +H03    Initialization
    !
    ! Sample Input
    !
    !  '+H03 TN408 Run PIN from inside DIVIMP  0-NO 1-YES        '    1
    !
    ! Input call documentation
    !
    !  'RUNPIN 0-NO 1-YES'
    !
    CPINOPT = 1
    !
    ! TAG:  +H04    Initialization
    !
    ! Sample Input
    !
    !  '+H04 ' 'TN408 Pin: reire07                              '
    !  READ(CPINCOM(11:80),'(A69)') ACTPIN
    !
    ! Input call documentation
    !
    !  'COMMAND TO RUN PIN'
    !
    CPINCOM = 'TN408 Pin: reire07                              '
    in = index(cpincom,':')
    READ(CPINCOM(in+1:),'(a69)') ACTPIN

    !
    ! TAG:  +H05    Initialization
    !
    ! Sample Input
    !
    !  '+H05      PIN Cell Area Option (IHCORR)                 '    1
    !
    ! Input call documentation
    !
    !  'PIN CELL AREA OPTION'
    !
    IHCORR = 1
    !
    ! TAG:  +H06    Initialization
    !
    ! Sample Input
    !
    !  '+H06      PIN Hybrid Wall Option 0=off 1,2=selection    '    0
    !
    ! Input call documentation
    !
    !  'HYBRID WALL IN PIN'
    !
    IHYBRID = 0
    !
    ! TAG:  +H07    Initialization
    !
    ! Sample Input
    !
    !  '+H07      PIN Puffing Option - 0=off 1=on               '    0
    !
    ! Input call documentation
    !
    !      'PIN PUFFING OPTION'
    !
    pinpuff = 0
    !
    ! TAG:  +H08    Initialization
    !
    ! Sample Input
    !
    !  '+H08      PIN Puff Location switch - 0=main SOL 1=PP    '    0
    !
    ! Input call documentation
    !
    !    'PUFF LOCATION OPTION'
    !
    swpvhpf = 0
    !
    ! TAG:  +H09    Initialization
    !
    ! Sample Input
    !
    !  '+H09      PIN Puff fraction (opt 1)                     '   0.0
    !
    ! Input call documentation
    !
    !  'PIN RE-PUFF FRAC'
    !
    hpcpuf = 0.0
    !
    ! TAG:  +H10    Initialization
    !
    ! Sample Input
    !
    !  '+H10      PIN Recycle -> Puff fraction (puff opt 2)     '   0.16
    !
    ! Input call documentation
    !
    !   'FLUX-> PUFF OPT 2'
    !
    ppcpuf = 0.16
    !
    ! TAG:  +H11    Initialization
    !
    ! Sample Input
    !
    !  '+H11      PIN Puff Injection temperature (eV)           '   0.5
    !
    ! Input call documentation
    !
    !  'PIN RE-PUFF TEMP '
    !
    tpufh = 0.5
    !
    ! TAG:  +H12    Initialization
    !
    ! Sample Input
    !
    !  '+H12      PIN Puff location indices JHPUF1(1 and 2)'  -17 -1000
    !
    ! Input call documentation
    !
    !  'PUFF LOCATION INDICES 1'
    !
    jhpuf1(1) = -17
    jhpuf1(2) = -1000
    !
    ! TAG:  +H13    Initialization
    !
    ! Sample Input
    !
    !  '+H13      PIN Puff location indices JHPUF2(1 and 2)'  -16 -1001
    !
    ! Input call documentation
    !
    !  'PUFF LOCATION INDICES 2'
    !
    jhpuf2(1) = -16
    jhpuf2(2) = -1001

    !
    ! TAG:  +H14    NIMBUS input block (A tag label is reserved for the NIMBUS block but tags are not used to read this input)
    !               A default input NAMELIST is assigned and is overwritten if one is read from the input file. 
    !
    ! NIMBUS input namelist
    !
    ! NOTE: This INPUT is not used by DIVIMP. If NIMBUS is being run it will access the input file and read this information directly using the
    !       namelist identifier and structure. Setting a default set and reading it in for tagged input files is mostly for documentation purposes
    !       in DIVIMP.
    !
    !       Set character array to the empty string
    
    cnimbin = ''

    ! ----------------------------------------------------------------------- 
    !
    !     HC related variable initializations 
    !
    !     TAG H15 to H64 and H90,H91: 
    !
    !     Insert code to initialize the HC variables 

    call global_hc_init 


    !  end subroutine initialize_tag_series_H


    !subroutine initialize_tag_series_I
    !implicit none

    !
    !  Initialization for TAG series I
    !
    !
    ! TAG:  +I01    Initialization
    !
    ! Sample Input
    !
    !  '+I01    Injection opt  1/2/3                            '     1
    !
    ! Input call documentation
    !
    !  'INJECTION OPT        '
    !
    CIOPTE = 1
    !
    ! TAG:  +I02    Initialization
    !
    ! Sample Input
    !
    !  '+I02    First diffuse  0inst 1random 2tpara             '     0
    !
    ! Input call documentation
    !
    !  'FIRST DIFFUSE OPT    '
    !
    CDIFOP = 0
    !
    ! TAG:  +I03    Initialization
    !
    ! Sample Input
    !
    !  '+I03    Control switch 0atoms 1ions                     '     0
    !
    ! Input call documentation
    !
    !  'CONTROL SWITCH       '
    !
    CNEUTA = 0
    !
    ! TAG:  +I04    Initialization
    !
    ! Sample Input
    !
    !  '+I04    Self- Sputtering Option 0-off 1-on              '     1
    !
    ! Input call documentation
    !
    !  'SELF-SPUTTER OPTION  '
    !
    CSELFS = 1
    !
    ! TAG:  +I05    Initialization
    !
    ! Sample Input
    !
    !  '+I05    Init ion Vel   1                                '     1
    !
    ! Input call documentation
    !
    !  'INITIAL ION VELOCITY '
    !
    CNEUTG = 1
    !
    ! TAG:  +I06    Initialization
    !
    ! Sample Input
    !
    !  '+I06 TN1465 Follow Imp.Ions Recombined to Neutrals 0=off '     1
    !
    ! Input call documentation
    !
    !  'FOLLOW REC. NEUTRAL  '
    !
    CFOLREC = 1
    !
    ! TAG:  +I07    Initialization
    !
    ! Sample Input
    !
    !  '+I07 TN1479 Ion Prompt Redeposition Option 0=off 1=on    '     0
    !
    ! Input call documentation
    !
    !  'ION PROMPT DEPOSITION'
    !
    prompt_depopt = 0
    !
    ! TAG:  +I08    Initialization
    !
    ! Sample Input
    !
    !  '+I08 T   Target Mirror Option 0-off  1-on                '     0
    !
    ! Input call documentation
    !
    !  'TARGET MIRROR OPT    '
    !
    cmiropt = 0
    !
    ! TAG:  +I09    Initialization
    !
    ! Sample Input
    !
    !  '+I09 T   Ion Periphery Option   0 to 3                   '     0
    !
    ! Input call documentation
    !
    !  'ION PERIPHERY OPT    '
    !
    FPOPT = 0
    !
    ! TAG:  +I10    Initialization
    !
    ! Sample Input
    !
    !  '+I10 TN996 Periphery Recycle Option       0-off 1-on     '     0
    !
    ! Input call documentation
    !
    !  'FP RECYCLE OPT       '
    !
    fpropt = 0
    !
    ! TAG:  +I11    Initialization
    !
    ! Sample Input
    !
    !  '+I11    Z effective (self)                 Zeff         '      1
    !
    ! Input call documentation
    !
    !  'ZEFF(SELF)     '
    !
    CIZEFF = 1
    !
    ! TAG:  +I12    Initialization
    !
    ! Sample Input
    !
    !  '+I12    Initial ionization state of impurity ions       '      1
    !
    ! Input call documentation
    !
    !  'INITIAL IZ STATE'
    !
    CIZSC = 1
    !
    ! TAG:  +I13    Initialization
    !
    ! Sample Input
    !
    !  '+I13    Collision Enhancement Factor       Zenh         '    1.0
    !
    ! Input call documentation
    !
    !  'COLLIS ENHANC ZENH'
    !
    CZENH = 1.0
    !
    ! TAG:  +I14    Initialization
    !
    ! Sample Input
    !
    !  '+I14    Set Ti = max(Ti,Tb) when reaching state (0 off) '      0
    !
    ! Input call documentation
    !
    !  'SET TI=TB AT STATE'
    !
    CIZSET = 0
    !
    ! TAG:  +I15    Initialization
    !
    ! Sample Input
    !
    !  '+I15    Maximum ionization state                        '      6
    !
    ! Input call documentation
    !
    !  'MAX IZ STATE'
    !
    NIZS = 6
    !
    ! TAG:  +I16    Initialization
    !
    ! Sample Input
    !
    !  '+I16    Stop following ions reaching Main Plasm 0no 1yes'      0
    !
    ! Input call documentation
    !
    !  'STOP WHEN HIT MAIN'
    !
    CSTOP  = 0
    !
    ! TAG:  +I17    Initialization
    !
    ! Sample Input
    !
    !  '+I17    Ion removal loss time              Tloss  (s)   '    0.000
    !
    ! Input call documentation
    !
    !  'ION LOSS TIME'
    !
    TLOSS = 0.000
    !
    ! TAG:  +I18    Initialization
    !
    ! Sample Input
    !
    !  '+I18 TN480 Ring for ion injection option 2        INJIR  '    1
    !
    ! Input call documentation
    !
    !  'RING NUMBER FOR INJ'
    !
    INJIR = 1
    !
    ! TAG:  +I19    Initialization
    !
    ! Sample Input
    !
    !  '+I19 TN480 Factor for Region Lower Bound          INJF1  '   0.0
    !
    ! Input call documentation
    !
    !  'INJECTION AREA LB' 
    !
    INJF1 = 0.0
    !
    ! TAG:  +I20    Initialization
    !
    ! Sample Input
    !
    !  '+I20 TN480 Factor for Region Upper Bound          INJF2  '   1.0
    !
    ! Input call documentation
    !
    !  'INJECTION AREA UB' 
    !
    INJF2 = 1.0
    !
    ! TAG:  +I21    Initialization
    !
    ! Sample Input
    !
    !  '+I21 TN443 X-max for Far-Periphery region (O/I) '       0.1  0.1
    !
    ! Input call documentation
    !
    !  'FP DIFFUSION SIZE' 
    !
    FPXMAXO = 0.1
    fpxmaxI = 0.1
    !
    ! TAG:  +I22    Initialization
    !
    ! Sample Input
    !
    !  '+I22 TN443 Far-Periphery Target Loss Time (O/I) '    1.0e-3  1.0e-3
    !
    ! Input call documentation
    !
    !  'FP TARGET LOSS T ' 
    !
    FPTIMO = 1.0e-3
    fptimI = 1.0e-3
    !
    ! TAG:  +I23    Initialization
    !
    ! Sample Input
    !
    !  '+I23 TN688 Far-Periphery Diffusion Rate ( < 0 = CDPERP ) '  -1.0
    !  IF (CDPERPFP.LT.0.0) CDPERPFP = CDPERP
    !
    ! Input call documentation
    !
    !  'FP DIFF RATE  ' 
    !
    CDPERPFP = -1.0

    ! ----------------------------------------------------------------------- 

    !     TAG I24 

    !     init_pos_opt - this option affects the initial position of neutrals 
    !     generated by ions which stike the target and subsequently recycle 
    !     as well as the initial position of ions that are formed from neutrals. 
    !     Option 1 for neutrals uses the cross component at the time of target 
    !     impact to estimate a more precise R,Z location for recycling - away 
    !     from the center of the target. For neutral ionization this option will 
    !     invoke getscross_approx to get a guesstimate of appropriate S CROSS 
    !     values for the initial position of the ion. 

    init_pos_opt = 1 

    ! ----------------------------------------------------------------------- 

    !     TAG I25 

    !     fp_neut_opt - this option turns on the possibility of neutral 
    !     ionization within a far peripheral region - it only applies to 
    !     wall elements not target ones. 
    !     = 0 = off 
    !     = 1 = on 

    fp_neut_opt = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG I26 

    !     fp_plasma_opt - this option specifies how a plasma will be 
    !     determined for the far periphery if a plasma is required. This 
    !     is particularly related to option I25. 
    !     = 0 = plasma from nearest associated real grid cell 
    !     = 1 = Te from specified input, ne from grid 
    !     = 2 = Both Te and Ne from specified input 

    fp_plasma_opt = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG I27 

    !     fp_te - default far periphery temperature if option in use (in eV) 

    fp_te = 10.0 

    ! ----------------------------------------------------------------------- 

    !     TAG I28 

    !     fp_ne - default far periphery density if option in use (m-3) 

    fp_ne = 1.0e18 

    ! ----------------------------------------------------------------------- 

    !     TAG I29 and I30 

    !     These are endpoints or corner points for the line/box ion injection 
    !     options ciopte=9 or ciopte=10. Default values are not 
    !     meaningful but are assigned to values of 0.0. Code for 
    !     ciopte 9 and 10 will check if all values are zero and will 
    !     issue and error and exit. 

    cxsca = 0.0 
    cysca = 0.0 
    cxscb = 0.0 
    cyscb = 0.0 

    ! ----------------------------------------------------------------------- 

    !     TAG I31 

    !     Far periphery transport flow option 
    !     0 - no flow in far periphery 
    !     1 - flow in far periphery is the same as associated ring 
    !     2 - flow in far periphery is specified as input 

    fp_flow_opt = 0 

    !     TAG I32 

    !     Far periphery variable to hold velocity input 

    fp_flow_velocity_input = 0.0 

    !------------------------------------------------------------------------ 

    !     TAG I33 

    !     Number of radial cells in simple far periphery mesh 

    fp_n_bins = 20 

    !     TAG I34 

    !     Defines the option used to choose the width of simple crude 
    !     FP mesh 
    !     Option 0 = maximum distance from edge cell to wall from 
    !               fp_walldist values 
    !     Option 1 = width of grid is specified using fpxmaxo for MAIN and 
    !                fpxmaxi for PFZ 

    fp_grid_width_opt = 0 



    !------------------------------------------------------------------------ 

    !     TAG I35 

    !     Ion reflection coefficients to roughly mimic ion pumping at the grid 
    !     edge in MC and PFZ 

    !     Main chamber - default is full reflection 

    mc_recyc = 1.0 

    !     TAG I36 

    !     Ion reflection coefficients to roughly mimic ion pumping at the grid 
    !     edge in MC and PFZ 

    !     Private Flux Zone - default is full reflection 

    pfz_recyc = 1.0 

    !------------------------------------------------------------------------ 

    !     TAG I37 

    !     Fluid code charge state index specifying the density distribution 
    !     to be used for injection option 12 and 13 

    !     NOTE: Depending on the fluid data in the file this value will not 
    !     necessarily be the same as the impurity charge state to be 
    !     followed. e.g. a SOLPS solution containing hydrogen, helium and 
    !     neon density data will load the hydrogen data into the non-impurity 
    !     arrays while the helium and neon are both loaded into the impurity 
    !     data arrays. If you want the ne+ density this is then index 3 since 
    !     there are data for he+ and he++ in the first two slots of the arrays. 

    !     Default value is 1 assuming you want to use the first charge state of 
    !     the first impurity as the basis for injection. 

    e2diz_inj = 1 

    ! TAG I38 
    ! Average charge state near the target for prompt deposition 
    ! option 4. The value is used in the gyroradius calculation to 
    ! determine if the ion promptly redeposits due to gryoradius 
    ! effects. 
    prompt_dep_avg_z = 1.0 

    !  end subroutine initialize_tag_series_I

    !      subroutine initialize_tag_series_K
    !      implicit none
    !------------------------------------------------------------------------ 

    !     TAG K?? 

    !     Series K unstructured input tags are assigned to the ERO-DIVIMP 
    !     interface code. These tags are both initialized and read by routines 
    !     contained in the ero_interface module 

    call init_ero_unstructured_input 

    !  end subroutine initialize_tag_series_K

    !    subroutine initialize_tag_series_N
    !      implicit none




    !
    !  Initialization for TAG series N
    !
    !
    ! TAG:  +N01    Initialization
    !
    ! Sample Input
    !
    !  '+N01    Launch option  0distrib 1point 2asymp 3tip 4wall'     0
    !
    ! Input call documentation
    !
    !  'LAUNCH OPTION        '
    !
    CNEUTB = 0
    !
    ! TAG:  +N02    Initialization
    !
    ! Sample Input
    !
    !  '+N02    Vel/angle flag 0-11                             '     0
    !
    ! Input call documentation
    !
    !  'VEL/ANGLE FLAG       '
    !
    CNEUTC = 0
    !
    ! TAG:  +N03    Initialization
    !
    ! Comments from input
    !
    !
    !     THESE VARIABLES (CNEUTH AND CNEUTI) WERE ADDED TO ALLOW LAUNCHES FROM THE WALLS,
    !     IN ADDITION TO THE NORMAL PLATE LAUNCHES. HOWEVER, THE
    !     CAPABILITY CAN BE USED TO SIMULATE TWO DIFFERING DISTRIBUTIONS
    !     OF PARTICLES COMING FROM THE PLATES. PERHAPS DUE TO
    !     PHYSICAL AND CHEMICAL SPUTTERING SIMULATNEOUSLY.
    !
    !
    ! Sample Input
    !
    !  '+N03 TN487 Supplemental Launch Option (as above)         '     0
    !
    ! Input call documentation
    !
    !  'SUP LAUNCH OPTION    '
    !
    CNEUTH = 0
    !
    ! TAG:  +N04    Initialization
    !
    ! Sample Input
    !
    !  '+N04 TN487 Supplemental V/A flag      (as above)         '     3
    !
    ! Input call documentation
    !
    !  'SUP VEL/ANGLE FLAG   '
    !
    CNEUTI = 3
    !
    ! TAG:  +N05    Initialization
    !
    ! Sample Input
    !
    !  '+N05    Initial Neutral Vel/Ang flag (-1=above,0-13)    '    -1
    !
    ! Input call documentation
    !
    !  'NEUT VEL/ANGLE FLAG  '
    !
    NVAOPT = -1
    !
    ! TAG:  +N06    Initialization
    !
    ! Sample Input
    !
    !  '+N06 TN1490 Supplemental 2D Neutral Launch 0=off 1=UEDGE '     0
    !
    ! Input call documentation
    !
    !  'EXTRA 2D LAUNCH OPTION'
    !
    neut2d_opt = 0
    !
    ! TAG:  +N07    Initialization
    !
    ! Sample Input
    !
    !  '+N07 TN1490 V/A Flag for 2D Neutral Launch (as regular)  '    3
    !
    ! Input call documentation
    !
    !  'EXTRA 2D V/A FLAG'
    !
    neut2d_vaopt = 3
    !
    ! TAG:  +N08    Initialization
    !
    ! Sample Input
    !
    !  '+N08    Sputter option 0std 1special 2mix 3self 4selfva1'     0
    !
    ! Input call documentation
    !
    !  'SPUTTER OPTION       '
    !
    CNEUTD = 0
    !
    ! TAG:  +N09    Initialization
    !
    ! Sample Input
    !
    !  '+N09 TN1209 Secondary Sputter Option                     '     0
    !
    ! Input call documentation
    !
    !  '2ND SPUTTER OPTION   '
    !
    CNEUTD2 = 0
    !
    ! TAG:  +N10    Initialization
    !
    ! Sample Input
    !
    !  '+N10    Normal option  0std 1fromT                      '     0
    !
    ! Input call documentation
    !
    !  'NORMAL OPTION        '
    !
    CNEUTE = 0
    !
    ! TAG:  +N11    Initialization
    !
    ! Sample Input
    !
    !  '+N11    NEUT spreading 0off 1on                         '     0
    !
    ! Input call documentation
    !
    !  'NEUT SPREADING       '
    !
    CNEUTF = 0
    !
    ! TAG:  +N12    Initialization
    !
    ! Sample Input
    !
    !  '+N12 TN1488 Neutral Impurity Velocity Type Option        '     0
    !
    ! Input call documentation
    !
    !  'IMP.NEUT.VEL.TYPE.OPT'
    !
    cneutvel = 0
    !
    ! TAG:  +N13    Initialization
    !
    ! Comments from input
    !
    ! ammod - added options 3 and 4 to NRFOPT
    !
    ! Sample Input
    !
    !  '+N13 T   Neutral Wall Reflection 0-off 1-on              '     0
    !
    ! Input call documentation
    !
    !  'NEUTRAL REFLECTION   '
    !
    NRFOPT = 0
    !
    ! TAG:  +N14    Initialization
    !
    ! Sample Input
    !
    !  '+N14 TN1481 Imp Neutral Momentum.Tran.Coll Opt 0=off 1=on'     2
    !
    ! Input call documentation
    !
    !  'NEUT.MOM.TRAN.COLL OPT '
    !
    mtcopt = 2
    !
    ! TAG:  +N15    Initialization
    !
    ! Sample Input
    !
    !  '+N15    Measure theta from T (degrees to X=0) for launch'   90.0
    !
    ! Input call documentation
    !
    !  'MEASURE THETA FROM'
    !
    CSNORM = 90.0
    !
    ! TAG:  +N16    Initialization
    !
    ! Sample Input
    !
    !  '+N16 ' 'TN487 Launch probability modifiers for each     '
    !  '    TN487 Wall segment range  #1  #2  mod.  :- '           0
    !
    ! Input call documentation
    !
    !  'WALL LAUNCH PROB. MODIFIERS'
    !
    !
    ! TAG:  +N17    Initialization
    !
    ! Sample Input
    !
    !  '+N17 TN721 Use Wall Probabilities as Absolute 0=No 1=Yes '    0
    !
    ! Input call documentation
    !
    !  'ABS WALL PROB    '
    !
    WLPABS = 0
    !
    ! TAG:  +N18    Initialization
    !
    ! Sample Input
    !
    !  '+N18 A18   Power of cosine release dist. (V/A 12,13)     '   1.0
    !
    ! Input call documentation
    !
    !  'COSINE DIST. POWER'
    !
    CNIN = 1.0
    !
    ! TAG:  +N19    Initialization
    !
    ! Sample Input
    !
    !  '+N19 TN???? Extra Velocity Multiplier for V/A 14         '   1.0
    !
    ! Input call documentation
    !
    !  'VELOCITY MULTIPLIER'
    !
    CVAMULT = 1.0
    !
    ! TAG:  +N20    Initialization
    !
    ! Sample Input
    !
    !  '+N20 TN???? Velocity Multiplier for Recombined Neutrals  '   1.0
    !
    ! Input call documentation
    !
    !  'REC. VEL. MULT'
    !
    CVRMULT = 1.0


    ! ----------------------------------------------------------------------- 

    !     TAG N21 : External sputtering flux data source 
    !               0 = geier file format for Ar 
    !               1 = import divimp charge resolved flux and energy 
    !                   data from a previous divimp run 

    !     Set the external flux data source to default to DIVIMP - add info in .dat file 
    !     in case someone wants to use the limited applicability Geier option 

    ext_flx_data_src = 1 

    !  end subroutine initialize_tag_series_N



    !    subroutine initialize_tag_series_P
    !      implicit none


    !
    !  Initialization for TAG series P
    !
    !
    ! TAG:  +P01    Initialization
    !
    ! Sample Input
    !
    !  '+P01    SOL option     0,1,1a,2,3,4,5,9,10  99file      '    22
    !
    ! Input call documentation
    !
    !  'SOL OPT              '
    !
    CIOPTF = 22
    !
    ! TAG:  +P02    Initialization
    !
    ! Sample Input
    !
    !  '+P02    Core Option    0,1,2,3                          '     0
    !
    ! Input call documentation
    !
    !  'Core Option        '
    !
    ccoreopt = 0
    !
    ! TAG:  +P03    Initialization
    !
    ! Sample Input
    !
    !  '+P03    Plasma decay   0std                 99file      '     0
    !
    ! Input call documentation
    !
    !  'PLASMA DECAY OPT     '
    !
    CIOPTG = 0
    !
    ! TAG:  +P04    Initialization
    !
    ! Comments from input
    !
    !
    !     Example
    !     10 23  3.0      4.0 22.0  0.0  0.0   0.0     3.0
    !     26 33  3.0      4.0 22.0  0.0  0.0   0.0     3.0
    !
    ! Sample Input
    !
    !  '+P04 ' 'BG PLASMA Options by Ring (PlasDec opts 90&91)  '
    !  '    R1,R2,Sect, PlasDec, Sol, Teg, Tig, Core, Efield'     0
    !
    ! Input call documentation
    !
    !  'SET OF BG PLASMA OPTIONS BY RING'
    !
    nbgplas=0

    !
    ! TAG:  +P05    Initialization
    !
    ! Sample Input
    !
    !  '+P05    Trap Tgrad Opt 0off 1on                         '     1
    !
    ! Input call documentation
    !
    !  'TRAP TGRAD OPTION    '
    !
    CIOPTO = 1
    !
    ! TAG:  +P06    Initialization
    !
    ! Sample Input
    !
    !  '+P06    SOL Enhancement Factor - Electric Field    SOLEF'    1.0
    !
    ! Input call documentation
    !
    !  'SOL  ENHANCE E(S) '
    !
    CSOLEF = 1.0
    !
    ! TAG:  +P07    Initialization
    !
    ! Sample Input
    !
    !  '+P07    SOL Enhancement Factor - Drift Velocity    SOLVF'    1.0
    !
    ! Input call documentation
    !
    !  'SOL  ENHANCE VH(S)'
    !
    CSOLVF = 1.0
    !
    ! TAG:  +P08    Initialization
    !
    ! Sample Input
    !
    !  '+P08    SOL1a Factor                               fl   '    0.01
    !
    ! Input call documentation
    !
    !  'SOL1A FACTOR "FL" '
    !
    CFL    = 0.01
    !
    ! TAG:  +P09    Initialization
    !
    ! Sample Input
    !
    !  '+P09    SOL1a Factor                               fs   '    1.0
    !
    ! Input call documentation
    !
    !  'SOL1A FACTOR "FS" '
    !
    CFS    = 1.0
    !
    ! TAG:  +P10    Initialization
    !
    ! Sample Input
    !
    !  '+P10    SOL10 Reversal Mach Number                 fRM  '    1.0
    !
    ! Input call documentation
    !
    !  'SOL10 FACTOR "FRM"'
    !
    CFRM   = 1.0
    !
    ! TAG:  +P11    Initialization
    !
    ! Sample Input
    !
    !  '+P11    SOL10 factor                               kin  '    1.0
    !
    ! Input call documentation
    !
    !  'SOL10 FACTOR "KIN"'
    !
    CKIN   = 1.0
    !
    ! TAG:  +P12    Initialization
    !
    ! Sample Input
    !
    !  '+P12    SOL10 factor                               kout '    1.2
    !
    ! Input call documentation
    !
    !  'SOL10 FACTOR"KOUT"'
    !
    CKOUT  = 1.2
    !
    ! TAG:  +P13    Initialization
    !
    ! Sample Input
    !
    !  '+P13    SOL10 factor                               fRmin'    0.01
    !
    ! Input call documentation
    !
    !  'SOL10 FACT "FRMIN"'
    !
    CFRMIN = 0.01
    !
    ! TAG:  +P14    Initialization
    !
    ! Sample Input
    !
    !  '+P14    SOL10 factor                               fRmax'    0.4
    !
    ! Input call documentation
    !
    !  'SOL10 FACT "FRMAX"'
    !
    CFRMAX = 0.4
    !
    ! TAG:  +P15    Initialization
    !
    ! Sample Input
    !
    !  '+P15    SOL6&7 Vb Length factor 1 (* SMAX)         VbL1 '    0.166
    !
    ! Input call documentation
    !
    !  'SOL11 LENGTH 1'
    !
    CVBL1  = 0.166
    !
    ! TAG:  +P16    Initialization
    !
    ! Sample Input
    !
    !  '+P16    SOL6&7 Vb multiplication factor 1          VbM1 '    0.0
    !
    ! Input call documentation
    !
    !  'SOL11 MULT 1  '
    !
    CVBM1  = 0.0
    !
    ! TAG:  +P17    Initialization
    !
    ! Sample Input
    !
    !  '+P17    SOL6&7 Vb Length factor 2 (* SMAX)         VbL2 '    0.5
    !
    ! Input call documentation
    !
    !  'SOL11 LENGTH 2'
    !
    CVBL2  = 0.5
    !
    ! TAG:  +P18    Initialization
    !
    ! Sample Input
    !
    !  '+P18    SOL6&7 Vb multiplication factor 2          VbM2 '    0.0
    !
    ! Input call documentation
    !
    !  'SOL11 MULT 2  '
    !
    CVBM2  = 0.0
    !
    ! TAG:  +P19    Initialization
    !
    ! Sample Input
    !
    !  '+P19    Power density                      P/A    (W/m2)'    3.0E+07
    !
    ! Input call documentation
    !
    !  'POWER DENSITY'
    !
    CPA   = 3.0E+07
    !
    ! TAG:  +P20    Initialization
    !
    ! Sample Input
    !
    !  '+P20    Parallel heat conduction coeff     K0           '    2.0E+03
    !
    ! Input call documentation
    !
    !  'PAR HEAT COND'
    !
    CK0   = 2.0E+03
    !
    ! TAG:  +P21    Initialization
    !
    ! Sample Input
    !
    !  '+P21    Parallel heat conduction -ions     K0I          '     58.9
    !
    ! Input call documentation
    !
    !  'PAR HEAT COND IONS'
    !
    CK0I  = 58.9
    !
    ! TAG:  +P22    Initialization
    !
    ! Sample Input
    !
    !  '+P22    Override input E-field from file E=0  0-off 1-on'      3
    !
    ! Input call documentation
    !
    !  'EFIELD=0 FOR FILE'
    !
    OFIELD = 3
    !
    ! TAG:  +P23    Initialization
    !
    ! Sample Input
    !
    !  '+P23 T1433 E-field Opt 4 - Length of E-field region *SMAX'    0.25
    !
    ! Input call documentation
    !
    !  'EFIELD SRC LENGTH'
    !
    CEFLEN = 0.25
    !
    ! TAG:  +P24    Initialization
    !
    ! Sample Input
    !
    !  '+P24 T1433 E-field Opt 4 - Te collisionality multiplier  '    2.0
    !
    ! Input call documentation
    !
    !  'EFIELD COLL-MULT'
    !
    CEFfact = 2.0
    !
    ! TAG:  +P25    Initialization
    !
    ! Sample Input
    !
    !  '+P25 TN401 Decay length for ionization source Ls  SOL12  '    0.08
    !
    ! Input call documentation
    !
    !   'LS DECAY SOL12   '
    !
    CSOLLS = 0.08
    !
    ! TAG:  +P26    Initialization
    !
    ! Sample Input
    !
    !  '+P26 T     Second decay characteristic length            '    0.08
    !
    ! Input call documentation
    !
    !   'SECOND DECAY DIST'
    !
    CSOLLT = 0.08
    !
    ! TAG:  +P27    Initialization
    !
    ! Sample Input
    !
    !  '+P27 T     Source fraction                               '    0.5
    !
    ! Input call documentation
    !
    !   'SOURCE FRACTION  '
    !
    CFIZ   = 0.5
    !
    ! TAG:  +P28    Initialization
    !
    ! Sample Input
    !
    !  '+P28 TN401 Decay length for radiative losses  Lr  SOL12  '    0.08
    !
    ! Input call documentation
    !
    !   'LR DECAY SOL12   '
    !
    CSOLLR = 0.08
    !
    ! TAG:  +P29    Initialization
    !
    ! Sample Input
    !
    !  '+P29 TN401 Coefficient for radiative losses   Pr  SOL12  '    1.0
    !
    ! Input call documentation
    !
    !   'PR CONST SOL12   '
    !
    CSOLPR = 1.0
    !
    ! TAG:  +P30    Initialization
    !
    ! Sample Input
    !
    !  '+P30 TN775 Radiation source strength fraction Frr        '    1.0
    !
    ! Input call documentation
    !
    !   'FR SOURCE FRAC.  '
    !
    CSOLFR = 1.0
    !
    ! TAG:  +P31    Initialization
    !
    ! Sample Input
    !
    !  '+P31 TN401 Source Ionization Option 0-lin 1-exp   SOL12  '    1
    !
    ! Input call documentation
    !
    !   'IONIZATION SOURCE OPT'
    !
    CSOPT  = 1
    !
    ! TAG:  +P32    Initialization
    !
    ! Sample Input
    !
    !  '+P32 TN401 Source Radiation Option  0-lin 1-exp   SOL12  '    3
    !
    ! Input call documentation
    !
    !   'RADIATION SOURCE OPT '
    !
    CPOPT  = 3
    !
    ! TAG:  +P33    Initialization
    !
    ! Sample Input
    !
    !  '+P33 TN660 Imaginary Root Option                  SOL12+ '    1
    !
    ! Input call documentation
    !
    !  'TREAT NEGATIVE ROOTS '
    !
    SROOTOPT = 1
    !
    ! TAG:  +P34    Initialization
    !
    ! Sample Input
    !
    !  '+P34 TN407 Flux Recirculation Option 0-off 1-on          '    0
    !
    ! Input call documentation
    !
    !  'SRC FLUX RECIRC OPT'
    !
    FLUXROPT = 0
    !
    ! TAG:  +P35    Initialization
    !
    ! Sample Input
    !
    !  '+P35 ' 'TN????+407 Set of flux specifications  '
    !  '    Ring , data1(I), data2(O), data3 : Rows - '          0
    !
    ! Input call documentation
    !
    !  'SET OF FLUX DATA'
    !
    FLUXPTS = 0

    !
    ! TAG:  +P36    Initialization
    !
    ! Sample Input
    !
    !  '+P36 TN408 Calculate SOL iteratively? 0-NO 1-YES         '    0
    !
    ! Input call documentation
    !
    !  'DO SOL AFTER PIN '
    !
    CITERSOL = 0
    !
    ! TAG:  +P37    Initialization
    !
    ! Sample Input
    !
    !  '+P37 TN408 Secondary SOL option for iterative calculation'   -2
    !  IF (CSECSOL.EQ.-2) CSECSOL = CIOPTF
    !
    ! Input call documentation
    !
    !  'SECONDARY SOL OPT'
    !
    CSECSOL  = -2
    !
    ! TAG:  +P38    Initialization
    !
    ! Sample Input
    !
    !  '+P38 TN408 Ionization option for iterative SOL           '    3
    !
    ! Input call documentation
    !
    !    'SECOND ION SOURCE OPT'
    !
    CSOPT2 = 3
    !
    ! TAG:  +P39    Initialization
    !
    ! Sample Input
    !
    !  '+P39      Number of Pin/SOL iterations                  '    3
    !
    ! Input call documentation
    !
    !  'NO. OF PIN ITER. '
    !
    NITERSOL = 3
    !
    ! TAG:  +P40    Initialization
    !
    ! Comments from input
    !
    !
    !     The following values apply to the specifiable SOL option that
    !     allows two part linearly interpolated fits for each of
    !     Te, Ti, ne and vb to be specified.
    !
    !
    ! The following lines specify the parameters for the linearly
    ! interpolated specified BG SOL option. This is purely empirical.
    !                    S-value      Function Value
    ! The form is:         0.0            F0
    !                       S1            F1
    !                       S2            F2
    ! For S>S2 F=F2
    !
    !
    ! Sample Input
    !
    !  '+P40 TN     Te S1 - First  S -value = ctes1 * SMAX       '    0.05
    !
    ! Input call documentation
    !
    !  'Te distance  1'
    !
    ctes1 = 0.05
    !
    ! TAG:  +P41    Initialization
    !
    ! Sample Input
    !
    !  '+P41 TN     Te F1 - First  Te-value = ctef1 * te0        '    1.5
    !
    ! Input call documentation
    !
    !  'Te fact/mult 1'
    !
    ctef1 = 1.5
    !
    ! TAG:  +P42    Initialization
    !
    ! Sample Input
    !
    !  '+P42 TN     Te S2 - Second S -value = ctes2 * SMAX       '    0.3
    !
    ! Input call documentation
    !
    !  'Te distance  2'
    !
    ctes2 = 0.3
    !
    ! TAG:  +P43    Initialization
    !
    ! Sample Input
    !
    !  '+P43 TN     Te F2 - Second Te-value = ctef2 * te0        '    2.0
    !
    ! Input call documentation
    !
    !  'Te fact/mult 2'
    !
    ctef2 = 2.0
    !
    ! TAG:  +P44    Initialization
    !
    ! Sample Input
    !
    !  '+P44 TN     Ti S1 - First  S -value = ctis1 * SMAX       '    0.05
    !
    ! Input call documentation
    !
    !  'Ti distance  1'
    !
    ctis1 = 0.05
    !
    ! TAG:  +P45    Initialization
    !
    ! Sample Input
    !
    !  '+P45 TN     Ti F1 - First  Te-value = ctif1 * ti0        '    1.5
    !
    ! Input call documentation
    !
    !  'Ti fact/mult 1'
    !
    ctif1 = 1.5
    !
    ! TAG:  +P46    Initialization
    !
    ! Sample Input
    !
    !  '+P46 TN     Ti S2 - Second S -value = ctis2 * SMAX       '    0.3
    !
    ! Input call documentation
    !
    !  'Ti distance  2'
    !
    ctis2 = 0.3
    !
    ! TAG:  +P47    Initialization
    !
    ! Sample Input
    !
    !  '+P47 TN     Ti F2 - Second Te-value = ctif2 * ti0        '    2.0
    !
    ! Input call documentation
    !
    !  'Ti fact/mult 2'
    !
    ctif2 = 2.0
    !
    ! TAG:  +P48    Initialization
    !
    ! Sample Input
    !
    !  '+P48 TN     Nb S1 - First  S -value = cnes1 * SMAX       '    0.05
    !
    ! Input call documentation
    !
    !  'Nb distance  1'
    !
    cnes1 = 0.05
    !
    ! TAG:  +P49    Initialization
    !
    ! Sample Input
    !
    !  '+P49 TN     Nb F1 - First  Te-value = cnef1 * ne0        '    2.0
    !
    ! Input call documentation
    !
    !  'Nb fact/mult 1'
    !
    cnef1 = 2.0
    !
    ! TAG:  +P50    Initialization
    !
    ! Sample Input
    !
    !  '+P50 TN     Nb S2 - Second S -value = cnes2 * SMAX       '    0.35
    !
    ! Input call documentation
    !
    !  'Nb distance  2'
    !
    cnes2 = 0.35
    !
    ! TAG:  +P51    Initialization
    !
    ! Sample Input
    !
    !  '+P51 TN     Nb F2 - Second Te-value = cnef2 * ne0        '    2.0
    !
    ! Input call documentation
    !
    !  'Nb fact/mult 2'
    !
    cnef2 = 2.0
    !
    ! TAG:  +P52    Initialization
    !
    ! Sample Input
    !
    !  '+P52 TN     vb S1 - First  S -value = cvbs1 * SMAX       '    0.1
    !
    ! Input call documentation
    !
    !  'vb distance  1'
    !
    cvbs1 = 0.1
    !
    ! TAG:  +P53    Initialization
    !
    ! Sample Input
    !
    !  '+P53 TN     vb F1 - First  Te-value = cvbf1 * ti0        '    0.25
    !
    ! Input call documentation
    !
    !  'vb fact/mult 1'
    !
    cvbf1 = 0.25
    !
    ! TAG:  +P54    Initialization
    !
    ! Sample Input
    !
    !  '+P54 TN     vb S2 - Second S -value = cvbs2 * SMAX       '    0.4
    !
    ! Input call documentation
    !
    !  'vb distance  2'
    !
    cvbs2 = 0.4
    !
    ! TAG:  +P55    Initialization
    !
    ! Sample Input
    !
    !  '+P55 TN     vb F2 - Second Te-value = cvbf2 * ti0        '    0.0
    !
    ! Input call documentation
    !
    !  'vb fact/mult 2'
    !
    cvbf2 = 0.0
    !
    ! TAG:  +P56    Initialization
    !
    ! Comments from input
    !
    !
    !     Marfe Core options - for option 4,5,6
    !
    !
    ! Sample Input
    !
    !  '+P56 TN1424 Core Option4,5- Velocity decay fraction *SMAX'    0.05
    !
    ! Input call documentation
    !
    !  'CORE-VEL FRAC'
    !
    corefv = 0.05
    !
    ! TAG:  +P57    Initialization
    !
    ! Sample Input
    !
    !  '+P57 TN1424 Core Option4,5- Te,Ti    decay fraction *SMAX'    0.05
    !
    ! Input call documentation
    !
    !  'CORE-TEMP FRAC'
    !
    coreft = 0.05
    !
    ! TAG:  +P58    Initialization
    !
    ! Sample Input
    !
    !  '+P58 TN1427 Core Option4,5- Velocity decay frac 2   *SMAX'    0.5
    !
    ! Input call documentation
    !
    !  'CORE-VEL FRAC2'
    !
    corefv2 = 0.5
    !
    ! TAG:  +P59    Initialization
    !
    ! Sample Input
    !
    !  '+P59 TN1427 Core Option4,5- Te,Ti    decay frac 2   *SMAX'    0.1
    !
    ! Input call documentation
    !
    !  'CORE-TEMP FRAC2'
    !
    coreft2 = 0.1






    ! ----------------------------------------------------------------------- 

    !     TAG P60 

    !     ngradopt - new density gradient option to allow for density 
    !     variation within the 2PM specification. Works in conjunction 
    !     with temperature gradient options. 

    ngradopt = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG P61 

    !     override_bg_velocity_opt - Option to override the background 
    !     velocity which is calculated by the other SOL options. 
    !     Option 0 : off 
    !     Option 1 : Prescribed flow - using data from osmns28 
    !     Option 2 : Recalculate background flow using the density 
    !                from the SOL option and source data from EIRENE 

    !     Default value is OFF 

    override_bg_velocity_opt = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG P62 

    !     ofield_targ - option for setting the target electric field value 
    !                   in the case of over-ride E-field calculations 

    !     option 1 : target E-field = 0 
    !                first cell data calculated as gradient to second cell only 
    !     option 2 : target E-field = first cell E-field 
    !                first cell data calculated as gradient to second cell only 
    !     option 3 : target E-field calculated from gradients to first cell center 
    !                first cell data calculated from gradients to target and to 
    !                second cell 

    !     default is option 3 

    ofield_targ = 3 


    ! ----------------------------------------------------------------------- 

    !     TAG P63 : External plasma overlay option 
    !               0 = off 
    !               1 = on 
    !               Default value is 0 

    external_plasma_overlay = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG P64 : External plasma overlay file name 
    !               - specifies the name of the file to be loaded 
    !               - full path required unless rundiv script is modified 
    !               Default value is null string 

    external_plasma_file = '' 

    ! ----------------------------------------------------------------------- 

    !     TAG P65 : SOL option 12/13 etc - additional pressure option 
    !               PADD - Adds additional pressure to (1+PADD) * PINF 
    !               Over a distance of PDIST * SMAX 

    sol13_padd = 0.0 

    !----------------------------------------------------------------------- 
    !     TAG P66 : SOL option 12/13 etc - additional pressure option 
    !               PDIST - Adds additional pressure to PMULT * PINF 
    !               Over a distance of PDIST * SMAX 

    sol13_pdist = 0.0 


    !  end subroutine initialize_tag_series_P




    !    subroutine intialize_tag_series_Q
    !      implicit none

    !
    !  Initialization for TAG series Q
    !
    !
    ! TAG:  +Q01    Initialization
    !
    ! Sample Input
    !
    !  '+Q01    TeB Gradient   0lin 1lin/lin 2p/a   99file      '     0
    !
    ! Input call documentation
    !
    !  'TEB GRADIENT OPTION  '
    !
    CIOPTK = 0
    !
    ! TAG:  +Q02    Initialization
    !
    ! Sample Input
    !
    !  '+Q02    TiB Gradient   0lin 1lin/lin 2p/a   99file      '     0
    !
    ! Input call documentation
    !
    !  'TIB GRADIENT OPTION  '
    !
    CIOPTL = 0
    !
    ! TAG:  +Q03    Initialization
    !
    ! Sample Input
    !
    !  '+Q03 TN???? Te,Ti Flatten Option                         '     0
    !
    ! Input call documentation
    !
    !  'Flatten Te
    !
    cflatopt = 0
    !
    ! TAG:  +Q04    Initialization
    !
    ! Sample Input
    !
    !  '+Q04 TN1447 Te flat upstream for S > SMAX * Teg cutoff = '    0.0
    !
    ! Input call documentation
    !
    !  'TEB GRAD. CUTOFF'
    !
    ctegcut = 0.0
    !
    ! TAG:  +Q05    Initialization
    !
    ! Sample Input
    !
    !  '+Q05 TN1447 Ti flat upstream for S > SMAX * Tig cutoff = '    0.0
    !
    ! Input call documentation
    !
    !  'TIB GRAD. CUTOFF'
    !
    ctigcut = 0.0
    !
    ! TAG:  +Q06    Initialization
    !
    ! Sample Input
    !
    !  '+Q06    Temperature of electrons at 0      TeB0  (eV)   '   20.0
    !
    ! Input call documentation
    !
    !  'TEMP AT 0   TEB0  '
    !
    CTEB0  = 20.0
    !
    ! TAG:  +Q07    Initialization
    !
    ! Sample Input
    !
    !  '+Q07    Temperature of electrons at plates TeBP  (eV)   '   20.0
    !
    ! Input call documentation
    !
    !  'PLATES TEMP TEBP  '
    !
    CTEBP  = 20.0
    !
    ! TAG:  +Q08    Initialization
    !
    ! Sample Input
    !
    !  '+Q08    Temperature outer TeB step         TeBout(eV)   '    0.0
    !
    ! Input call documentation
    !
    !  'TEMP STEP   TEBOUT'
    !
    CTEBOU = 0.0
    !
    ! TAG:  +Q09    Initialization
    !
    ! Sample Input
    !
    !  '+Q09    Temperature inner TeB step         TeBin (eV)   '  100.0
    !
    ! Input call documentation
    !
    !  'TEMP STEP   TEBIN '
    !
    CTEBIN = 100.0
    !
    ! TAG:  +Q10    Initialization
    !
    ! Sample Input
    !
    !  '+Q10    Temperature of trapped plasma      TeBt  (eV)   '   10.0
    !
    ! Input call documentation
    !
    !  'TEMP TRAP   TEBT  '
    !
    CTEBT  = 10.0
    !
    ! TAG:  +Q11    Initialization
    !
    ! Sample Input
    !
    !  '+Q11 TN1278 Te exp decay step in trap      TeBouP       '    0.01
    !
    ! Input call documentation
    !
    !  'TEMP DECAY TEBOUTP'
    !
    CTEBOUP = 0.01
    !
    ! TAG:  +Q12    Initialization
    !
    ! Sample Input
    !
    !  '+Q12    Temperature gradient factor        feBL1        '    0.0
    !
    ! Input call documentation
    !
    !  'GRAD FACTOR FEBL1 '
    !
    CFEBL1 = 0.0
    !
    ! TAG:  +Q13    Initialization
    !
    ! Sample Input
    !
    !  '+Q13    Temperature gradient factor        feBL2        '    0.0
    !
    ! Input call documentation
    !
    !  'GRAD FACTOR FEBL2 '
    !
    CFEBL2 = 0.0
    !
    ! TAG:  +Q14    Initialization
    !
    ! Sample Input
    !
    !  '+Q14    Temperature gradient factor        feBt         '    1.0
    !
    ! Input call documentation
    !
    !  'GRAD FACTOR FEBT  '
    !
    CFEBT  = 1.0
    !
    ! TAG:  +Q15    Initialization
    !
    ! Sample Input
    !
    !  '+Q15    Temperature gradient factor        feB2         '    1.0
    !
    ! Input call documentation
    !
    !  'GRAD FACTOR FEB2  '
    !
    CFEB2  = 1.0
    !
    ! TAG:  +Q16    Initialization
    !
    ! Sample Input
    !
    !  '+Q16    Temperature of ions at 0           TiB0  (eV)   '   20.0
    !
    ! Input call documentation
    !
    !  'TEMP AT 0   TIB0  '
    !
    CTIB0  = 20.0
    !
    ! TAG:  +Q17    Initialization
    !
    ! Sample Input
    !
    !  '+Q17    Temperature of ions at plates      TiBP  (eV)   '   20.0
    !
    ! Input call documentation
    !
    !  'PLATES TEMP TIBP  '
    !
    CTIBP  = 20.0
    !
    ! TAG:  +Q18    Initialization
    !
    ! Sample Input
    !
    !  '+Q18    Temperature outer TiB step         TiBout(eV)   '    0.05
    !
    ! Input call documentation
    !
    !  'TEMP STEP   TIBOUT'
    !
    CTIBOU = 0.05
    !
    ! TAG:  +Q19    Initialization
    !
    ! Sample Input
    !
    !  '+Q19    Temperature inner TiB step         TiBin (eV)   '  250.0
    !
    ! Input call documentation
    !
    !  'TEMP STEP   TIBIN '
    !
    CTIBIN = 250.0
    !
    ! TAG:  +Q20    Initialization
    !
    ! Sample Input
    !
    !  '+Q20    Temperature of trapped plasma      TiBt  (eV)   '   10.0
    !
    ! Input call documentation
    !
    !  'TEMP TRAP   TIBT  '
    !
    CTIBT  = 10.0
    !
    ! TAG:  +Q21    Initialization
    !
    ! Sample Input
    !
    !  '+Q21 TN1278 Ti exp decay step in trap      TiBouP       '    0.01
    !
    ! Input call documentation
    !
    !  'TEMP DECAY TIBOUTP'
    !
    CTIBOUP = 0.01
    !
    ! TAG:  +Q22    Initialization
    !
    ! Sample Input
    !
    !  '+Q22    Temperature gradient factor        fiBL1        '    0.0
    !
    ! Input call documentation
    !
    !  'GRAD FACTOR FIBL1 '
    !
    CFIBL1 = 0.0
    !
    ! TAG:  +Q23    Initialization
    !
    ! Sample Input
    !
    !  '+Q23    Temperature gradient factor        fiBL2        '    0.0
    !
    ! Input call documentation
    !
    !  'GRAD FACTOR FIBL2 '
    !
    CFIBL2 = 0.0
    !
    ! TAG:  +Q24    Initialization
    !
    ! Sample Input
    !
    !  '+Q24    Temperature gradient factor        fiBt         '    1.0
    !
    ! Input call documentation
    !
    !  'GRAD FACTOR FIBT  '
    !
    CFIBT  = 1.0
    !
    ! TAG:  +Q25    Initialization
    !
    ! Sample Input
    !
    !  '+Q25    Temperature gradient factor        fiB2         '    1.0
    !
    ! Input call documentation
    !
    !  'GRAD FACTOR FIB2  '
    !
    CFIB2  = 1.0
    !
    ! TAG:  +Q26    Initialization
    !
    ! Sample Input
    !
    !  '+Q26    Density at 0                       NB0   (m**-3)'    1.0E19
    !
    ! Input call documentation
    !
    !  'DENS AT 0   NB0   '
    !
    CNB0   = 1.0E19
    !
    ! TAG:  +Q27    Initialization
    !
    ! Sample Input
    !
    !  '+Q27 TN1264 Density at plates              NEBP  (m**-3)'    1.0e19
    !
    ! Input call documentation
    !
    !  'PLATES DENS NEBP  '
    !
    CNEBP  = 1.0e19
    !
    ! TAG:  +Q28    Initialization
    !
    ! Sample Input
    !
    !  '+Q28    Density outer NB step              NBout (m**-3)'    0.03
    !
    ! Input call documentation
    !
    !  'DENS STEP   NBOUT '
    !
    CNBOUT = 0.03
    !
    ! TAG:  +Q29    Initialization
    !
    ! Sample Input
    !
    !  '+Q29    Density inner NB step              NBin  (m**-3)'    5.0E18
    !
    ! Input call documentation
    !
    !  'DENS STEP   NBIN  '
    !
    CNBIN  = 5.0E18
    !
    ! TAG:  +Q30    Initialization
    !
    ! Sample Input
    !
    !  '+Q30    Density of trapped plasma          NBt   (m**-3)'    1.0E19
    !
    ! Input call documentation
    !
    !  'DENS TRAP   NBT   '
    !
    CNBT   = 1.0E19
    !
    ! TAG:  +Q31    Initialization
    !
    ! Sample Input
    !
    !  '+Q31 TN1278 Ne exp decay step in trap      NBouP        '    0.01
    !
    ! Input call documentation
    !
    !  'DENS DECAY  CNBOUP'
    !
    CNBOUP = 0.01
    !
    ! TAG:  +Q32    Initialization
    !
    ! Comments from input
    !
    !     ENTER LANGMUIR PROBE DATA ALONG PLATES FOR TEMPERATURE GRADIENT
    !     OPTIONS 3 AND 4, AND PLASMA DECAY OPTIONS 2 AND 3. IN THE
    !     CASE OF THE PLASMA DECAY OPTIONS THE DATA ENTERED HERE WILL
    !     BE FILLED IN FOR ALL POINTS OF THE PLASMA. THUS GIVING A
    !     UNIFORM CONDITION UNLESS A TEMPERATURE GRADIENT OPTION
    !     MODIFIES THE TEMPERATURE.
    !
    !     FOR THOSE OPTIONS WITH ONLY ONE SET OF DATA THE INNER
    !     SPECIFICATION IS USED.
    !
    !     NOTE LPDATI :- LANGMUIR PROBE DATA INNER ...
    !
    !     THE LANGMUIR PROBE DATA SWITCH IDENTIFIES THE THIRD QUANTITY
    !     AS EITHER THE DENSITY OR THE Isat (PROBE SATURATION CURRENT)
    !
    !     THE DATA IS ORDERED AS  RING #,TEBP,TIBP,NBP (or ISAT)
    !
    !
    ! Sample Input
    !
    !  '+Q32 TN1347 Langmuir Probe Switch     0=Nb  1=Isat       '     1
    !
    ! Input call documentation
    !
    !  'LP DATA SWITCH'
    !
    lpdatsw = 1
    !
    ! TAG:  +Q33    Initialization
    !
    ! Sample Input
    !
    !  '+Q33    OUTER Target Data Multipliers (Te,Ti,Nb):'  1.0   1.0   1.0
    !
    ! Input call documentation
    !
    !  'Inner Multipliers       '
    !
    te_mult_i = 1.0
    ti_mult_i = 1.0
    n_mult_i = 1.0
    !
    ! TAG:  +Q34    Initialization
    !
    ! Comments from input
    !
    !
    !     Outer
    !
    !
    ! Sample Input
    !
    !  '+Q34 ' 'Probe data at outer plate (opt4) or total (opt3)'
    !  '    Ring , TeBP , TiBP , NBP : Number of rows - '          0
    !
    ! Input call documentation
    !
    !  'SET OF L. PROBE DATA INNER'
    !

    NLPDATI = 0


    !
    ! TAG:  +Q35    Initialization
    !
    ! Sample Input
    !
    !  '+Q35    INNER Target Data Multipliers (Te,Ti,Nb):'  1.0   1.0   1.0
    !
    ! Input call documentation
    !
    !  'Outer Multipliers       '
    !
    te_mult_o = 1.0
    ti_mult_o = 1.0
    n_mult_o = 1.0
    !
    ! TAG:  +Q36    Initialization
    !
    ! Sample Input
    !
    !  '+Q36 ' 'Probe data at inner plate          (T grad opt4)'
    !  '    Ring , TeBP , TiBP , NBP : Number of rows - '          0
    !
    ! Input call documentation
    !
    !  'SET OF L. PROBE DATA OUTER'
    !
    NLPDATO = 0


    !
    ! TAG:  +Q37    Initialization
    !
    ! Comments from input
    !
    !
    !     Read in data for core rings - data varies depending on
    !     core option selected.
    !
    !
    ! Sample Input
    !
    !  '+Q37 ' 'CORE Plasma Data - for Core Options 1,2 and 3'
    !  '    Ring , TeB , TiB , NB , Vb : Number of rows - '       0
    !
    ! Input call documentation
    !
    !  'SET OF CORE DATA'
    !
    Ncoredat = 0

    !
    ! TAG:  +Q38    Initialization
    !
    ! Sample Input
    !
    !  '+Q38    Inboard plasma flow velocity       Vhin  (m/s)  '    0.0
    !
    ! Input call documentation
    !
    !  'INB. PLASMA FLOW V'
    !
    CVHIN = 0.0
    !
    ! TAG:  +Q39    Initialization
    !
    ! Sample Input
    !
    !  '+Q39    Inboard electric field             Ein   (V/m)  '    0.0
    !
    ! Input call documentation
    !
    !  'INBOARD ELEC FIELD'
    !
    CEIN  = 0.0
    !
    ! TAG:  +Q40    Initialization
    !
    ! Sample Input
    !
    !  '+Q40    Outboard plasma flow vel  (SOL5)   Vhout (m/s)  '    0.0
    !
    ! Input call documentation
    !
    !  'OUT. PLASMA FLOW V'
    !
    CVHOUT = 0.0
    !
    ! TAG:  +Q41    Initialization
    !
    ! Sample Input
    !
    !  '+Q41    Outboard electric field   (SOL5)   Eout  (V/m)  '    0.0
    !
    ! Input call documentation
    !
    !  'OUTBRD ELEC FIELD '
    !
    CEOUT  = 0.0


    ! ----------------------------------------------------------------------- 

    !     TAG Q42 

    !     Tags Q43 and Q43 specify a temperature on a ring by ring 
    !     basis which is then used instead of the target temperature for 
    !     calculating the target heat flux and sputtering yields. At the 
    !     present time these values are used directly in NEUT. 

    !     TAG Q42 

    !     - Two parameters - IR TE - for INNER JET/OUTER SONNET 

    nsheath_vali = 0 
    !      call rzero(sheath_vali,maxnrs*2) 


    ! ----------------------------------------------------------------------- 

    !     TAG Q43 

    !     - Two parameters - IR TE - for OUTER JET/INNER SONNET 

    nsheath_valo = 0 
    !      call rzero(sheath_valo,maxnrs*2) 


    ! ----------------------------------------------------------------------- 

    !     TAG Q44 - Core plasma profiles as a function of PSIN 

    !     - Five parameters - PSIN TE TI NE VB 

    !     Data array is allocatable and does not need initialization 
    !     (note: use allocatable_data module 

    ncoreprofile = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG Q45 - Shift applied to psin coordinate of core profiles 

    delta_psin_core = 0.0 

    !  end subroutine intialize_tag_series_Q


    !    subroutine initialize_tag_series_R
    !      implicit none

    !
    !  Initialization for TAG series R
    !
    !
    ! TAG:  +R01    Initialization
    !
    ! Sample Input
    !
    !  '+R01      DETACHED PLASMA: Length Scaling Switch 0-S 1-P'    0
    !
    ! Input call documentation
    !
    !  'S21 Length Ref Switch'
    !
    s21refsw = 0
    !
    ! TAG:  +R02    Initialization
    !
    ! Sample Input
    !
    !  '+R02 TN988 Detached Plasma Model: Te/Te0 at L1 O/I '  1.0    1.0
    !
    ! Input call documentation
    !
    !  'Te Ratio at L1 O/I'
    !
    terat = 1.0
    terati = 1.0
    !
    ! TAG:  +R03    Initialization
    !
    ! Sample Input
    !
    !  '+R03 TN1439   Ti/Ti0 at Position L1            O/I '  1.0    1.0
    !
    ! Input call documentation
    !
    !  'Ti Ratio at L1 O/I'
    !
    tirat = 1.0
    tirati = 1.0
    !
    ! TAG:  +R04    Initialization
    !
    ! Sample Input
    !
    !  '+R04 TN988    Ne/Ne0 at Position L1            O/I ' 10.0   10.0
    !
    ! Input call documentation
    !
    !  'Ne Ratio at L1 O/I'
    !
    nrat = 10.0
    nrati = 10.0
    !
    ! TAG:  +R05    Initialization
    !
    ! Sample Input
    !
    !  '+R05 TN1496   Exponent for Ne Equation 1.0=lin O/I '  1.0    1.0
    !
    ! Input call documentation
    !
    !  'Ne Exponent at L1 O/I'
    !
    nalph = 1.0
    nalphi = 1.0
    !
    ! TAG:  +R06    Initialization
    !
    ! Sample Input
    !
    !  '+R06 TN988    Qrad/Q0                          O/I '  5.0    5.0
    !
    ! Input call documentation
    !
    !  'Emitted energy rat'
    !
    qrat = 5.0
    qrati = 5.0
    !
    ! TAG:  +R07    Initialization
    !
    ! Sample Input
    !
    !  '+R07 TN988    L1/SMAX ratio                    O/I '  0.1    0.1
    !
    ! Input call documentation
    !
    !  'L1*SMAX boundary 1'
    !
    l1rat = 0.1
    l1rati = 0.1
    !
    ! TAG:  +R08    Initialization
    !
    ! Sample Input
    !
    !  '+R08 TN988    L2/SMAX ratio                    O/I '  0.2    0.2
    !
    ! Input call documentation
    !
    !  'L2*SMAX boundary 2'
    !
    l2rat = 0.2
    l2rati = 0.2
    !
    ! TAG:  +R09    Initialization
    !
    ! Sample Input
    !
    !  '+R09 TN988    LV/SMAX ratio                    O/I '  0.2    0.2
    !
    ! Input call documentation
    !
    !  'V->0 at lvrat*SMAX'
    !
    lvrat = 0.2
    lvrati = 0.2
    !
    ! TAG:  +R10    Initialization
    !
    ! Sample Input
    !
    !  '+R10 TN       Velocity multiplier at L1        O/I '  1.0    1.0
    !
    ! Input call documentation
    !
    !  'Vmult for Region A'
    !
    vbmult = 1.0
    vbmulti = 1.0
    !
    ! TAG:  +R11    Initialization
    !
    ! Sample Input
    !
    !  '+R11  ' 'TN    DETACHED Plasma Specifications on a by ring basis'
    !  'TN    INNER - IR TER TIR NR QR L1R L2R LVR VBM - N= '    0
    !
    ! Input call documentation
    !
    !  'SOL21 PARAMS INNER'
    !
    ns21i=0


    !
    ! TAG:  +R12    Initialization
    !
    ! Sample Input
    !
    !  '+R12 ' 'TN    DETACHED Plasma Specifications on a by ring basis'
    !  'TN    OUTER - IR TER TIR NR QR L1R L2R LVR VBM - N= '    0
    !
    ! Input call documentation
    !
    !  'SOL21 PARAMS OUTER'
    !
    ns21o=0


    ! ----------------------------------------------------------------------- 

    !     TAG R13 

    !     The following R-tags are enhancements to the detached plasma 
    !     prescription that allow for the shape of the density profile 
    !     in this region to be more precisely specified. This was done 
    !     to make it possible to specify profiles in this region which 
    !     correspond reasonably well to the measured Divertor Thomson 
    !     profiles. By specifying values for the following inputs - ONLY 
    !     supported on a ring by ring basis - the default standard 
    !     behaviour is used elsewhere - the density profiles in region A 
    !     of the detached plasma prescription can be more generally 
    !     defined. 

    !     R13 - array of input data for INNER JET/OUTER Sonnet 

    !     TAG R13 

    !     - 9 parameters - IR L1A L1B NR1A NR1B TER1A TER1B TIR1A TIR1B 

    aux_ns21i = 0 
    !      call rzero(aux_s21parmi,maxnrs*9) 


    ! ----------------------------------------------------------------------- 

    !     TAG R14 

    !     R14 - array of input data for OUTER JET/INNER Sonnet 

    !     - 9 parameters - IR L1A L1B NR1A NR1B TER1A TER1B TIR1A TIR1B 

    aux_ns21o = 0 
    !      call rzero(aux_s21parmo,maxnrs*9) 

    ! ----------------------------------------------------------------------- 

    !     TAG R15 

    !     R15 - aux velocity factor 


    aux_vel21 = 1.0 

    !  end subroutine initialize_tag_series_R

    !    subroutine initialize_tag_series_S
    !      implicit none


    !
    !  Initialization for TAG series S
    !
    !
    ! TAG:  +S01    Initialization
    !
    ! Sample Input
    !
    !  '+S01    On-AXIS Toroidal B-field value                  '    2.0
    !
    ! Input call documentation
    !
    !  'ON-AXIS B-FIELD'
    !
    cbphi = 2.0
    !
    ! TAG:  +S02    Initialization
    !
    ! Sample Input
    !
    !  '+S02    Mass of plasma ions                Mb           '    2.0
    !
    ! Input call documentation
    !
    !  'PLASMA ION MASS  '
    !
    CRMB = 2.0
    !
    ! TAG:  +S03    Initialization
    !
    ! Sample Input
    !
    !  '+S03    Charge on plasma ions              Zb           '    1
    !  RIZB = REAL (CIZB)
    !
    ! Input call documentation
    !
    !  'PLASMA ION CHARGE'
    !
    CIZB = 1
    RIZB = REAL (CIZB)
    !
    ! TAG:  +S04    Initialization
    !
    ! Sample Input
    !
    !  '+S04    Mass of impurity ions              Mi           '   12.0
    !
    ! Input call documentation
    !
    !  'IMPURITY ION MASS'
    !
    CRMI = 12.0
    !
    ! TAG:  +S05    Initialization
    !
    ! Sample Input
    !
    !  '+S05    Atomic number of impurity ions     Zi           '      6
    !
    ! Input call documentation
    !
    !  'IMP ATOMIC NUMBER'
    !
    CION = 6
    !
    ! TAG:  +S06    Initialization
    !
    ! Sample Input
    !
    !  '+S06    Initial temperature                Tem1  (eV)   '    0.5
    !
    ! Input call documentation
    !
    !  'INITIAL TEMP   '
    !
    CTEM1 = 0.5
    !
    ! TAG:  +S07    Initialization
    !
    ! Sample Input
    !
    !  '+S07    Initial temperature (2)            Tem2  (eV)   '    0.0
    !
    ! Input call documentation
    !
    !  'INITIAL TEMP (2)'
    !
    CTEM2  = 0.0
    !
    ! TAG:  +S08    Initialization
    !
    ! Sample Input
    !
    !  '+S08    Initial R position of impurity     R0    (m)    '    1.4
    !
    ! Input call documentation
    !
    !  'INITIAL R POSITION'
    !
    CXSC = 1.4
    !
    ! TAG:  +S09    Initialization
    !
    ! Sample Input
    !
    !  '+S09    Initial Z position of impurity     Z0    (m)    '    0.9
    !
    ! Input call documentation
    !
    !  'INITIAL Z POSITION'
    !
    CYSC = 0.9
    !
    ! TAG:  +S10    Initialization
    !
    ! Comments from input
    !
    !     jdemod - remove upper bounds checks for maxizs and maximp due
    !              to dynamical allocation
    !
    ! Sample Input
    !
    !  '+S10    Operation Mode  1 Time-Dependent  2 Steady-State'      2
    !
    ! Input call documentation
    !
    !  'OPERATION MODE'
    !
    IMODE = 2
    !
    ! TAG:  +S11    Initialization
    !
    ! Sample Input
    !
    !  '+S11    Number of impurity ions to be followed          '    100
    !
    ! Input call documentation
    !
    !  'NO OF IONS'
    !
    NIMPS = 100
    !
    ! TAG:  +S12    Initialization
    !
    ! Sample Input
    !
    !  '+S12 TN487 Number of Supplementary Neutrals to Launch   '      0
    !
    ! Input call documentation
    !
    !  'NUM SUP IONS'
    !
    NIMPS2 = 0
    !
    ! TAG:  +S13    Initialization
    !
    ! Sample Input
    !
    !  '+S13    Quantum iteration time for atoms   fsrate (s)   '    1.0E-8
    !
    ! Input call documentation
    !
    !  'NEUT ITERATE TIME'
    !
    FSRATE = 1.0E-8
    !
    ! TAG:  +S14    Initialization
    !
    ! Sample Input
    !
    !  '+S14    Quantum iteration time for ions    qtim   (s)   '    1.0E-8
    !
    ! Input call documentation
    !
    !  'DIV ITERATE TIME'
    !
    QTIM = 1.0E-8
    !
    ! TAG:  +S15    Initialization
    !
    ! Comments from input
    !
    !---- READ IN TIME DEPENDENT DATA  (NOTE 128)
    !
    !
    ! Sample Input
    !
    !  '+S15 T   CPU time limit in seconds          cpulim (s)   '  80000.0
    !
    ! Input call documentation
    !
    !  'CPU TIME LIMIT'  
    !
    CPULIM = 80000.0
    !
    ! TAG:  +S16    Initialization
    !
    ! Sample Input
    !
    !  '+S16 ' 'Average Dwell Times (s) for each state 0,1,2..'
    !  '        Number of dwell times given below :-'    0
    !
    ! Input call documentation
    !
    !  'DWELT'
    !
    NQS = 0
    !
    ! TAG:  +S17    Initialization
    !
    ! Sample Input
    !
    !  '+S17 ' 'Dwell Time Factors for time-dependent analysis'
    !  '        Number of dwell time factors :-'  0
    !
    ! Input call documentation
    !
    !  'T FACTORS'
    !
    NTS = 0
    !
    ! TAG:  +S18    Initialization
    !
    ! Sample Input
    !
    !  '+S18 T   Maximum dwell time for steady state (s)         '    0.5
    !
    ! Input call documentation
    !
    !  'DWELL TIME LIMIT '
    !
    ctimmax = 0.5
    !
    ! TAG:  +S19    Initialization
    !
    ! Sample Input
    !
    !  '+S19    Random number seed  (0 generate new seed)       '      0
    !
    ! Input call documentation
    !
    !  'RANDOM NUMBER SEED'
    !
    CISEED = 0
    !
    ! TAG:  +S20    Initialization
    !
    ! Sample Input
    !
    !  '+S20    Number of Iterations                            '      1
    !
    ! Input call documentation
    !
    !  'NO. OF ITERATIONS '
    !
    NITERS = 1
    !
    ! TAG:  +S21    Initialization
    !
    ! Sample Input
    !
    !  '+S21    SOLTEST - 0.0 run normally -1.0 test SOL opt    '     0.0
    !
    ! Input call documentation
    !
    !  'TEST SOL OPT ONLY'
    !
    ctestsol = 0.0





    ! ----------------------------------------------------------------------- 

    !     TAG S22 

    !     S22 - Option to turn on neutral V/A flag debugging 
    !           OFF by default 

    debug_neutv = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG S23 

    !     S23 - maximum energy/velocity used in the neutral velocity 
    !           debugging options. 

    debug_neutv_einmax = 50.0 

    ! ----------------------------------------------------------------------- 

    !     TAG S24 

    !     S24 - Number of bins to divide the velocity distribution when 
    !           debugging neutral velocity 

    debug_neutv_nbins = 500 


    !  end subroutine initialize_tag_series_S

    !    subroutine intialize_tag_series_T
    !    implicit none

    !
    !  Initialization for TAG series T
    !
    !
    ! TAG:  +T01    Initialization
    !
    ! Sample Input
    !
    !  '+T01    Ionization     0/1/2old 3ei 4no 5ei/dis 6no/dis '     3
    !
    ! Input call documentation
    !
    !  'IONIZATION OPT       '
    !
    CIOPTA = 3
    !
    ! TAG:  +T02    Initialization
    !
    ! Comments from input
    !
    !     All other Reiser input is unstructured and optional
    !
    ! Sample Input
    !
    !  '+T02    Collision opt  0std 1inf 2zeff 3tpara           '    13
    !
    ! Input call documentation
    !
    !  'COLLISION OPT        '
    !
    CIOPTB = 13
    !
    ! TAG:  +T03    Initialization
    !
    ! Sample Input
    !
    !  '+T03    REISER                                          '     0
    !
    ! Input call documentation
    !
    !  'REISER COLLISION OPT '
    !
    CIOPTR = 0
    !
    ! TAG:  +T04    Initialization
    !
    ! Sample Input
    !
    !  '+T04    Friction opt   0std 1inf 2tpara                 '     0
    !
    ! Input call documentation
    !
    !  'FRICTION OPT         '
    !
    CIOPTC = 0
    !
    ! TAG:  +T05    Initialization
    !
    ! Sample Input
    !
    !  '+T05    Heating opt    0std 1inf 2zero 3Ti              '     0
    !
    ! Input call documentation
    !
    !  'HEATING OPT          '
    !
    CIOPTD = 0
    !
    ! TAG:  +T06    Initialization
    !
    ! Sample Input
    !
    !  '+T06    CX Recomb opt  0off 1on 2Vcx                    '     0
    !
    ! Input call documentation
    !
    !  'CX RECOMB OPT        '
    !
    CIOPTI = 0
    !
    ! TAG:  +T07    Initialization
    !
    ! Sample Input
    !
    !  '+T07    Dperp option   0const 1vary                     '     0
    !
    ! Input call documentation
    !
    !  'DPERP OPTION         '
    !
    CIOPTJ = 0
    !
    ! TAG:  +T08    Initialization
    !
    ! Sample Input
    !
    !  '+T08 TN1272 Perpendicular step option 0-normal 1-core    '     3
    !
    ! Input call documentation
    !
    !  'PERP STEP PROB OPTION'
    !
    cdiffopt = 3
    !
    ! TAG:  +T09    Initialization
    !
    ! Sample Input
    !
    !  '+T09 TN14?? Pinch Velocity Option 0=off 1=all 2=main SOL '     0
    !
    ! Input call documentation
    !
    !  'PINCH VELOCITY OPTION'
    !
    pinchopt = 0
    !
    ! TAG:  +T10    Initialization
    !
    ! Sample Input
    !
    !  '+T10    TeB Grad Coeff 0off 1on                         '     1
    !
    ! Input call documentation
    !
    !  'TEB GRAD COEFF OPTION'
    !
    CIOPTM = 1
    !
    ! TAG:  +T11    Initialization
    !
    ! Sample Input
    !
    !  '+T11    TiB Grad Coeff 0off 1on                         '     1
    !
    ! Input call documentation
    !
    !  'TIB GRAD COEFF OPTION'
    !
    CIOPTN = 1
    !
    ! TAG:  +T12    Initialization
    !
    ! Sample Input
    !
    !  '+T12    T-GRAD Force Modification Function 0-off 1-UEDGE'     0
    !
    ! Input call documentation
    !
    !  'GRAD FORCE MOD OPTION'
    !
    fgradopt = 0
    !
    ! TAG:  +T13    Initialization
    !
    ! Sample Input
    !
    !  '+T13 TN505 Poloidal Velocity Drift Option 0-off 1-on     '     0
    !
    ! Input call documentation
    !
    !  'POL. DRIFT OPTION    '
    !
    CPDRFT = 0
    !
    ! TAG:  +T14    Initialization
    !
    ! Sample Input
    !
    !  '+T14    Cross Field Diffusion factor       Dperp (m*m/s)'    0.3
    !
    ! Input call documentation
    !
    !  'X DIFF DPERP       '
    !
    CDPERP = 0.3
    !
    ! TAG:  +T15    Initialization
    !
    ! Comments from input
    !
    ! move after input file read
    !
    ! Sample Input
    !
    !  '+T15    Trap Cross Field Diffusion factor  Dperpt(m*m/s)'   -1.0
    !  if (cdperpt.le.0.0) cdperpt = cdperp
    !
    ! Input call documentation
    !
    !  'X DIFF DPERP-TRAP'
    !
    CDPERPT = -1.0
    !
    ! TAG:  +T16    Initialization
    !
    ! Sample Input
    !
    !  '+T16    Perpendicular Pinch Velocity        Vpinch (m/s)'    0.0
    !
    ! Input call documentation
    !
    !  'PERP PINCH VEL.  '
    !
    cVPINCH = 0.0
    !
    ! TAG:  +T17    Initialization
    !
    ! Sample Input
    !
    !  '+T17 TN505 Poloidal Drift Velocity (m/s)                 '   0.0
    !
    ! Input call documentation
    !
    !  'POL. DRIFT VEL.  '
    !
    CDRFTV = 0.0
    !
    ! TAG:  +T18    Initialization
    !
    ! Sample Input
    !
    !  '+T18 TN    Poloidal Drift Range F1*SMAX < S < F2*SMAX'    0.2 0.8
    !
    ! Input call documentation
    !
    !  'Start and End of DRFTV'
    !
    cdrftv_start = 0.2
    cdrftv_end = 0.8


    ! ----------------------------------------------------------------------- 

    !     TAG T19 TO T27 

    !     Initialization of Reiser options to default values 

    aswitch = 0 
    sk11 = 0 
    sk12 = 0 
    sk13 = 0 
    sd11 = 0 
    sd12 = 0 
    sd13 = 0 
    coulomb_log = 15.0 
    linearpeak = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG T28 

    !     T28 - pinch_loc_opt - this option specifies the region of the 
    !                           grid where the radial velocity in 
    !                           pinchopt 4 will be applied. 

    !         = 0 = Entire grid excluding PFZ (default) 
    !         = 1 = main SOL only 
    !         = 2 = Entire grid including PFZ 
    !               - this requires a sign change to the value assigned to 
    !                 Vr (or Vpinch depending on terminology) 
    !         = 3 = Main SOL only above Xpoint region 
    !               (i.e. In cells adjacent to adjacent to Xpoint and above) 
    !         = 4 = Main SOL above Xpoint region + core 


    pinch_loc_opt = 0 


    ! ----------------------------------------------------------------------- 

    !     TAG T29 

    !     T29 - pinch_npdf, pinch_pdf - this option loads the probability 
    !           distribution function to be used when randomly 
    !           determining the value of the pinch/radial velocity at 
    !           each time step. 

    !           At the present time no default PDF is loaded - a check must 
    !           be added at the end of the read routine to make sure 
    !           that a PDF has been specified if pinchopt 4 has been selected. 

    pinch_npdf = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG T30 

    !     T30 - pinch correlation time. A new radial velocity value will be 
    !           chosen periodically based on the value of this quantity. 
    !           The default value of 0.0 will result in a new velocity 
    !           being selected every time step. This value is specified 
    !           in second. On JET it is typically 5 to 20 microseconds. 

    pinch_correlation_time = 0.0 

    ! ----------------------------------------------------------------------- 

    !     TAG T31 

    !     T31 - Drift region - specifies the region to which poloidal 
    !           drifts should be applied. 
    !           1 - SOL + PFZ 
    !           2 - SOL only 
    !           3 - PFZ only 
    !           4 - CORE only 

    !           Other options can easily be added as needed - the default is 
    !           option 1. 

    drft_region = 1 

    ! ----------------------------------------------------------------------- 

    !     TAG T32 

    !     T32 - Drift Mach Option - Detailed drift velocity input on a ring 
    !           ring basis is specified as a mach number to be multiplied 
    !           by the sound speed at the top of the torus for each ring 

    !           Option 0 : OFF 
    !           Option 1 : CS calculated from 2*Te 
    !           Option 2 : CS calculated from Te+Ti 

    !           Default value is 0 - OFF - data is specified in terms of 
    !                                      velocity 

    drftvel_machopt = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG T33 

    !     T33 - Detailed specifications of data ring by ring - this is 
    !           an array listing 
    !                      ring number       velocity/mach 
    !           Data does not need to be specified for each ring - the 
    !           default value will be applied instead. 

    !           The number of array elements is initialized to zero 

    ndrftvel = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG T34 

    !     T34 - S displacement in 2D resulting from a perpendicular step in 
    !           paramagnetic "Z" direction in 3D - this actually moves the 
    !           particle onto an adjacent flux tube - however, since 
    !           DIVIMP is 2D - the effect is to actually move the particle 
    !           onto an adjacent identical flux tube at a different S location 
    !           - thus effectively giving a net S displacement. 

    !           The first approximation to this is to use 

    !           ds = cross_step * Btor/Bpol 

    !           Where Btot = sqrt(Btor**2+Bpol**2) 

    !           Bpol = sqrt(Br**2 + Bz**2) 

    !           A value of 0 for this option is OFF 
    !                      1 is ON 

    !           Default is shown below 

    dperpz_opt = 0 
    base_dperpz_step = 0.0 

    ! ----------------------------------------------------------------------- 

    !     TAG T35 - related to poloidal drift options - T31,T32,T33 

    !     T35 - Drift region calculation option 

    !     Option 0: Input values are specified in terms of S 
    !     Option 1: Input values are specified in terms of P (poloidal distance) 
    !     Option 2: Input values are specified in terms of Z (single null only) 

    !     Default value is S para 

    drft_distopt = 0 

    ! ----------------------------------------------------------------------- 

    !     TAG T36 to T39 - related to radial and poloidal ExB drift options 

    !     Default values for options related to calculating plasma potential 
    !     and exb drifts 

    !     TAG T36 
    !             - This option is used to determine the method of calculatng 
    !               plasma potential 
    !             - potopt = 0    Use 3xTe(0) at each target as the floating 
    !                              potential start start point 
    !             - potopt = 1    Import LP data listing the measured floating 
    !                              potential ... if imported data not available it 
    !                              defaults to option 0.(not yet implemented) 

    potopt = 0 

    !     TAG T37 
    !     exb_rad_opt = 0 ... no exb radial drift is applied 
    !                 = 1 ... exb radial drift is turned on 

    exb_rad_opt = 0 

    !     TAG T38 
    !     exb_pol_opt = 0 ... no exb poloidal drift is applied 
    !                 = 1 ... exb poloidal drift is turned on 

    exb_pol_opt = 0 

    !     TAG 39 

    !     exb_scale - real number 
    !               - the basic function of this is to switch the sign of 
    !                 the ExB drift for cases of forward (+1.0) and reverse (-1.0) 
    !                 B-field orientation. However, it can also be used as a scaling 
    !                 factor if the drifts are found to be either too large or too 
    !                 small. 

    exb_scale = 1.0 


    ! ----------------------------------------------------------------------- 

    !     Force scaling factors - all default to 1.0 
    !     T40 to T44 
    !     T40 = Friction force scaling factor (SF_FRIC) 
    !     T41 = Ion temperature force scaling factor (SF_TI) 
    !     T42 = Electron temperature force scaling factor (SF_TE) 
    !     T43 = Electric field force scaling factor (SF_EF) 
    !     T44 = Velocity diffusion scaling factor (SF_VDIFF) 
    !     T45 = Scaling factor for TAU (sf_tau) 

    !     Defaults: 
    !     sf_fric = 1.0 
    !     sf_ti   = 1.0 
    !     sf_te   = 1.0 
    !     sf_ef   = 1.0 
    !     sf_vdiff= 1.0 
    !     sf_tau  = 1.0 

    !----------------------------------------------------------------------- 
    sf_fric = 1.0 
    sf_ti   = 1.0 
    sf_te   = 1.0 
    sf_ef   = 1.0 
    sf_vdiff= 1.0 
    sf_tau  = 1.0 

    ! ----------------------------------------------------------------------- 

    !    T46 Velocity based temperature calculation option 

    !    0  = default 
    !    1+ = other options 

    ti_calc_opt = 0 


    ! ----------------------------------------------------------------------- 

    !    T47 Coulomb logarithm calculation options 

    !     0  = default = constant (default value = 15.0) 
    !     1  = 30.0 - 0.5 * LOG(ni) + 1.5 * LOG(ti)  [HC code - ] 
    !          Originally in Sivukhin, D.V., Coulomb collisions in a fully ionized plasma in 
    !          Review of Plasma Physics (Consultation Bureau, New York, 1966) Vol. 4, p.88. 

    !     2  = 17.3 - 0.5*LOG(n/1.0E20) + 1.5*LOG(t/1000.0)  [LIM code] 
    !     3  = log(1.5e13 * t**(1.5) / sqrt(n))   [SOL22 PEI term] 

    !     Default constant value 

    lambda_opt = 0 

    !    T48 Coulomb logarithm calculation options 
    !     Coulomb logarithm constant value - default value is 15.0 - this allows 
    !     specification of alternate constant values for option 0. 

    !     dafault value = 15.0 

    lambda_val = 15.0 

    ! T49 Blob frequency. Determines whether or not a radial pinch 
    ! velocity is selected. See description in unstructered_input_com. 
    ! The default value of -1 sets it to 1/qtim in div.f. That is, 
    ! a blob is encountered every time step. 
    fblob = -1 

    ! T50 Core diffusion coefficient. If -1.0 then cdperpc = cdperp. 
    cdperpc = -1.0 

    ! T51-53 are free to use. They were old options I scrapped. 

    ! T54: Approximate ballooning transport by the factor BT/BT @ OMP. 
    balloon_opt = 0 

    ! T55: For VR-PDF option, multiply chosen radial velocity values 
    ! by this when in the divertor. 
    div_vr_fact = 1.0 

    ! T56: Turn off parallel transport when impurity within blob 
    ! (determined by pinch_correlation_time). 
    in_blob = .false. 
    in_blob_switch = 0 

    ! T57: Turn on hole-like transport. The model assumes holes and 
    ! blobs are birthed at a given R-Rsep @ OMP (input via T61) each 
    ! with the same frequency at that location. As one moves outwards, 
    ! the frequency of holes exponentially decays according to a given 
    ! 1/e falloff length (T58). 
    hole_switch = 0.0 

    ! T58: The exponential decay length for the hole frequency as one 
    ! moves inwards in meters. Positive = decays as one moves OUTWARDS. 
    ! Not used if T57 = 0. 
    hole_lambda = 1.0 

    ! T59: Additional inward pinch velocity for the the core region 
    ! only. This operates independently of the other pinch options, so 
    ! no matter what those assign this pinch velocity is added on 
    ! after the fact in the core region only. 
    core_pinch = 0.0 

    ! T60: The minimum R-Rsep @ OMP value (in meters) for which the 
    ! blob/hole-like impurity transport model is used. The default is 
    ! 0.0 (SOL only) but one can set this to anything, e.g., -0.005 
    ! to allow blob/hole-like transport into the core a bit. 
    blob_min_rmrsomp = 0.0 

    ! T61: Birth location of hole/blob pairs as R-Rsep @ OMP (m). 
    ! At this location the frequency of blobs and holes are both 
    ! "fblob", the value input for T49. The physical meaning is that 
    ! for every blob, a hole is likewise born. 
    ! Outwards of this value: 
    !   fblob(r) = fblob 
    !   fhole(r) = fblob * exp(-r/hole_lambda) 
    ! Inwards of this value: 
    !   fblob(r) = fblob * exp(r/blob_lambda) 
    !   fhole(r) = fblob 
    ! Where 'r' is R-Rsep @ OMP. If hole-like transport (T57) is not 
    ! on then fhole(r) = 0.0 always. Defaults to 0.0, i.e., the 
    ! separatrix. 
    blob_birth_rmrsomp = 0.0 

    ! T62: The exponential decay length for the blob frequency as one 
    ! moves inwards in meters. Positive = decays as one moves INWARDS. 
    hole_lambda = 1.0 


    !  end subroutine initialize_tag_series_T



    !    subroutine initialize_tag_series_W
    !      implicit none
    ! ----------------------------------------------------------------------- 

    !     TAG W01 

    !     wall_plasma_opt - this option specifies the algorithm to be used 
    !     to define the plasma conditions (if any) associated with each 
    !     element of the wall. The target elements are associated with the 
    !     target plasma conditions and are not affected by this option. 

    !     0 = use plasma conditions in associated cell 
    !     1 = linear decay beyond last ring + interpolation along ring 
    !     2 = exponential decay beyond last ring + interpolation along ring 

    wall_plasma_opt = 0 


    ! ----------------------------------------------------------------------- 

    !     TAG W02 

    !     wall_plasma_fact - this is a scale factor to be used in the scaling 
    !     algorithms defined by wall_plasma_opt = 1,2 - it is set to 0.1 m for 
    !     now. There are minimum values for wall plasma conditions specified 
    !     in the code. These could be moved to optional input if required in 
    !     the future. 

    wall_plasma_fact = 0.1 

    ! ----------------------------------------------------------------------- 

    !  end subroutine initialize_tag_series_W

    return 
  end subroutine initializeunstructuredinput



  ! ====================================================================== 

  ! subroutine: ValidateUnstructuredInput 


  SUBROUTINE ValidateUnstructuredInput 
    use mod_params 
    use mod_comtor 
    IMPLICIT none 
    !     include 'params' 
    !     include 'comtor' 

    !integer :: maxchk 
    !PARAMETER (MAXCHK = 15) 
    !integer :: i1,i2,nchklist,mchklist(maxchk) 
    !character*72 :: chklist(2*maxchk) 
    character*69 :: sp 


    ! jdemod - default values have been assigned for all of these inputs so none should be required at this point
    !        - in addition, the checking code appears to not be run at all anyway since it is at an unreachable line label
    !
    !...  List of required unstructured input tags: 
    !DATA chklist &
    !     /'026',' PIN selection  0-NIMBUS 1-EIRENE97 2-EIRENE99     (1)', &
    !     '010',' Geometry data  0-standard 1-from DIVIMP              ', &
    !     '021',' Input file     0-standard 1-from DIVIMP              ', &
    !     '020','   Run time (CPU seconds)                             ', &
    !     '022','   Material: target  1-Mo 2-C 3-W 4-Be                ', &
    !     '024','             wall                                     ', &
    !     '011','   Grid type     0-structured 1-generalized           ', &
    !     '018','   Wall data     0-standard   1-seamless              ', &
    !     '019','   Debug option  0-off                                ', &
    !     'E11',' n-n collisions  0-off 1-standard mesh             (0)', &
    !     
    !     'E12',' Lyman alpha opacity  0-off 1-rec 2-rec&ion        (0)', &
    !     '058',' 1.1 Pressure gauge specification:                    ', &
    !     '076',' 1.0 Surface properties:                              ', &
    !     '077',' 1.0 Additional surfaces:                             ', &
    !     '078','   Target data shift  JET-i/o CMOD,DIIID-o/i (m) (0.0)'/ 


    ! jdemod 

    !     For unstructured input items whose default value is intended to match a 
    !     regular DIVIMP input - assign the DIVIMP values if the unstructured input 
    !     still contains its default value. At present this applies to some of the HC 
    !     code input data. 

    call global_hc_assign_inputs 

    ! jdemod 


    ! jdemod 

    !     For ion injection options 9 and 10 make sure that values have been 
    !     specified for the endpoints/corners of the injection line/region 

    if (cneuta == 1.and.(ciopte == 9.or.ciopte == 10).and. &
         (cxsca == 0.0.and.cysca == 0.0.and. &
         cxscb == 0.0.and.cyscb == 0.0)) then 

       !         Invalid inputs specified 

       write(0,*) 'ERROR: Invalid Input - INJECTION OPTION =', ciopte 
       write(0,*) '       INJECTION REGION NOT SPECIFIED' 
       write(0,*) '       SEE INPUT ITEMS *I29 and *I30' 
       write(6,*) 'ERROR: Invalid Input - INJECTION OPTION =', ciopte 
       write(6,*) '       INJECTION REGION NOT SPECIFIED' 
       write(6,*) '       SEE INPUT ITEMS *I29 and *I30' 
       write(7,*) 'ERROR: Invalid Input - INJECTION OPTION =', ciopte 
       write(7,*) '       INJECTION REGION NOT SPECIFIED' 
       write(7,*) '       SEE INPUT ITEMS *I29 and *I30' 

    end if

    if (cneuta == 0.and.(cneutb == 6.or.cneutb == 7).and. &
         (cxsca == 0.0.and.cysca == 0.0.and. &
         cxscb == 0.0.and.cyscb == 0.0)) then 

       !         Invalid inputs specified 

       write(0,*) 'ERROR: Invalid Input- NEUT LAUNCH OPTION =',cneutb 
       write(0,*) '       NEUTRAL LAUNCH REGION NOT SPECIFIED' 
       write(0,*) '       SEE INPUT ITEMS *I29 and *I30' 

       write(6,*) 'ERROR: Invalid Input- NEUT LAUNCH OPTION =',cneutb 
       write(6,*) '       NEUTRAL LAUNCH REGION NOT SPECIFIED' 
       write(6,*) '       SEE INPUT ITEMS *I29 and *I30' 

       write(7,*) 'ERROR: Invalid Input- NEUT LAUNCH OPTION =',cneutb 
       write(7,*) '       NEUTRAL LAUNCH REGION NOT SPECIFIED' 
       write(7,*) '       SEE INPUT ITEMS *I29 and *I30' 

    end if

    ! jdemod 

    RETURN 
    !
    ! jdemod - apparently unused input checking code below this point - inaccessible anyway so commenting it out
    !

!90  FORMAT(A,A3,A) 
!97  DO i1 = 1, MAXCHK 
!       if (mchklist(i1) /= 1) then 
!          WRITE(6,90) '  ',chklist(2*i1-1),chklist(2*i1) 
!          WRITE(0,90) '  ',chklist(2*i1-1),chklist(2*i1) 
!       end if
!    end do
!    GOTO 99 
!98  WRITE(6,90) 'SAMPLE: ',chklist(2*i2-1),chklist(2*i2) 
!99  STOP 
  END subroutine validateunstructuredinput




  !
  !==============================================================================
  !==============================================================================
  ! READ Routines
  !==============================================================================
  !==============================================================================
  !
  !
  !     The following routine reads in the unstructured input (aka optional 
  !     input) in the OEDGE/DIVIMP input file. The optional inputs are identified 
  !     by a tag in the input file. These tags are a letter followed by a 2-digit 
  !     number - leaving room for up to 2600 options in the current configuration 
  !     These options are usually arranged into letter groupings that are related 
  !     to a common theme. In addition, "tag series" which have a large number of 
  !     options are read in by tag specific subroutines: 

  !     e.g. 

  !        SUBROUTINE ReadTagSeries_G 
  !        SUBROUTINE ReadTagSeries_H 
  !        SUBROUTINE ReadTagSeries_I 

  !     In addition to the letter indexed tags the original implementation used 
  !     three digit numbers as the tags. These options are somewhat less well 
  !     documented. 

  !     Finally, the OEDGE/DIVIMP input file has a tag assigned for every input 
  !     value whether optional or not. This indexing will allow for all of the input 
  !     to be made optional when this is desirable. Any additional optional input 
  !     must be assigned a new and unique tag. 


  SUBROUTINE ReadUnstructuredInput(line2) 
    USE mod_osm_input 
    use allocatable_input_data 
    use mod_sol22_input 
    use ero_interface 
    use mod_params 
    use mod_slcom 
    use mod_cadas 
    use mod_comtor 
    use mod_cgeom 
    use mod_cedge2d 
    use mod_solparams 
    use mod_solswitch 
    !use mod_solcommon 
    use mod_reiser_com 
    use mod_line_profile 
    !use mod_out_unstruc 
    use mod_driftvel 
    use mod_diagvel 
    use mod_dperpz 
    use mod_lambda 
    use mod_sol23_input
    use mod_sol29_input 
    use mod_io 
    use mod_rundiv_local
    use mod_pindata
    use mod_dynam4
    use mod_adpak_com
    use mod_reader
    IMPLICIT none 


    !     jdemod - cleaning up some of the line length constraints since the inputs really don't need them 

    !      CHARACTER line2*(*),LINE*72,TAG*3,coment*72,cdum1*1024 
    character :: line2*(*),line*128,tag*3,cdum1*1024 
    real :: r,vol,z1,version 
    integer :: i,ir,ierr,i1,i2 

    !real :: getinputversion 

    ! jdemod - added variable to hold line read when calling RDG1 to get 
    !          ADAS data. 

    character :: line3*512 

    integer :: tagnum,fp 
    DATA    tagnum /0/ 

    COMMON /OPTTEMP/ osm_matcht,forcet1 
    integer :: osm_matcht,forcet1 

    !COMMON /MACHCOM/ machine2 
    !character*64 :: machine2 

    integer :: idum1 
    real :: rdum1 

    integer :: itag 
    logical :: status,sol28_first 
    character :: local_buffer*1024 

    DATA sol28_first /.TRUE./ 
    SAVE 

    !     Function declaration for TAG T29 

    !      real vtest,res,vr_pdf_int 
    !      external vr_pdf_int 
    integer :: in 
    real :: deltav1, deltav2 

    WRITE(line,'(A128)') line2 

    WRITE(TAG,'(A3)') LINE(3:5) 
    ! jdemod - add support for case insensitive tags
    call upcase(tag)

    ierr = 0 

    fp = PINOUT 

    ntaglist = ntaglist + 1 
    taglist(ntaglist) = tag 


    WRITE(SLOUT,*) 'TAG:',tag 


    !write(0,*) 'RUI START:',trim(tag),':',trim(line),':'
    

    !     *** NOTE *** 

    !...  This routine is getting unruly so it is being divided into 
    !     a series of routines that process the various input categories 
    !     individually: 
    ! ----------------------------------------------------------------------- 
    if     (line(2:2) == '{') then 
       IF (sol28_first) THEN 
          !          WRITE(0,*) 'INITIALIZING SOL28 OPTIONS' 
          !          CALL InitializeOptions 
          sol28_first = .FALSE. 
       end if
       !...    SOL28/OSM input options: 
       if (line(3:5) == '999'.or.line(3:6) == 'exit') then   ! not great... 
          CALL ProcessIterationBlocks 
       ELSE 
          s28mode = 4.0 
          WRITE(local_buffer,'(A)') line2 
          !          WRITE(0,*) 'local_buffer:',local_buffer(1:100) 
          !...      Isolate tag string: 
          DO itag = 2, LEN_TRIM(local_buffer) 
             if (local_buffer(itag:itag) == '}') exit 
          end do
          !...      Remove the portion of the data line that comes after a comment 
          !         character: 
          DO i1 = 1, LEN_TRIM(local_buffer) 
             if (local_buffer(i1:i1) == '$'.or. &
                  local_buffer(i1:i1) == '*') exit 
          end do
          local_buffer(i1:LEN_TRIM(local_buffer)) = ' ' 
          CALL ProcessInputTag(5,itag,local_buffer,status) 
       end if
       ntaglist = ntaglist - 1 

    else if (tag(1:1) == '2') then 
       ! ----------------------------------------------------------------------- 

       !     TAG 2??: 200 series tags are related to SOL22. 

       !     200 series options are unstructured input for SOL22 - the routine to process 
       !         these is in mod_sol22_input 

       call sol22_unstructured_input(tag,line,ierr) 


    else if (tag(1:1) == '3') then 
       ! ----------------------------------------------------------------------- 

       !     TAG 2??: 300 series tags are related to SOL23. 

       !     300 series options are unstructured input for SOL23 - the routine to process 
       !         these is in mod_sol23_input 

       call sol23_read_unstructured_input(tag,line,ierr) 


    else if (tag(1:1) == 'G') then 
       ! Series G 
       CALL ReadTagSeries_G(tag,line,fp) 


    else if (tag(1:1) == 'H') then 
       ! Series H 
       CALL ReadTagSeries_H(line,tag,fp) 


    else if (tag(1:1) == 'I') then 
       ! Series I 
       CALL ReadTagSeries_I(line,tag,fp) 


    else if (tag(1:1) == 'K') then 

       ! Series K 
       !     - this tag series is assigned to ERO interface related 
       !       quantities 

       CALL read_ero_unstructured_input(line,tag,fp) 

    else if (tag(1:1) == 'X') then 
       ! TAG X: Options related to SOL29. 
       call sol29_unstructured_input(tag, line) 

       ! ----------------------------------------------------------------------- 






       !
       ! Tag series 0 is Steve's original code and the tags are apparently in a semi-random order
       !

       ! Read series 0  (zero) unstructured input
       !
       ! Steve's intial unstructured input code - may need some modification to work with ypdates
       ! ----------------------------------------------------------------------- 
       !      ELSEIF (TAG(1:3).EQ.'001') THEN 
       ! ----------------------------------------------------------------------- 
       !      ELSEIF (TAG(1:3).EQ.'002') THEN 
       ! ----------------------------------------------------------------------- 
       !      ELSEIF (TAG(1:3).EQ.'003') THEN 
       ! ----------------------------------------------------------------------- 
       !      ELSEIF (TAG(1:3).EQ.'004') THEN 
       ! ----------------------------------------------------------------------- 
       !      ELSEIF (TAG(1:3).EQ.'005') THEN 
       ! ----------------------------------------------------------------------- 
       !      ELSEIF (TAG(1:3).EQ.'006') THEN 
       ! ----------------------------------------------------------------------- 
       !      ELSEIF (TAG(1:3).EQ.'007') THEN 
       ! ----------------------------------------------------------------------- 
    elseif (tag(1:3) == '008') then 
       CALL ReadI(line,outmode,0,3,'Output mode for user information') 
       ! ----------------------------------------------------------------------- 
    else if (tag(1:3) == '009') then 
       CALL ReadI(line,cgridst,0,1,'Stuff grid') 
       ! ----------------------------------------------------------------------- 

       !     EIRENE related options: 

    else if (tag(1:3) == '010') then 
       CALL ReadI(line,eirgeom,0,3,'EIRENE geometry') 
    else if (tag(1:3) == '011') then 
       CALL ReadI(line,eirgrid,0,1,'EIRENE grid type') 
    else if (tag(1:3) == '014') then 
       CALL ReadI(line,eiradd,0,100,'EIRENE AddUsr option') 
    else if (tag(1:3) == '018') then 
       CALL ReadI(line,eirneut,0,1,'EIRENE neutral data option') 
    else if (tag(1:3) == '019') then 
       CALL ReadI(line,eirdebug,-100,10000,'EIRENE debugging option') 
    else if (tag(1:3) == '020') then 
       CALL ReadI(line,eirtime,0,600000,'EIRENE execution time') 

       CALL GetEnv('DIVNAME',machine2) 
       !        WRITE(0,*) 'MARK: MACHINE2= '//machine2(1:LEN_TRIM(machine2)) 

       !       IPP/01 - Krieger: SUN Fortran chokes if machine has zero chars 

       if (.false..and.len_trim(machine2) > 0) then 
          if     (machine2(1:len_trim(machine2)) == 'hannah') then 
             CALL WN('GetInput','Increasing EIRENE runtime for HANNAH') 
             eirtime = eirtime * 1.1 
          else if (machine2(1:len_trim(machine2)) == 'jonah') then 
             !            CALL WN('GetInput','Decreasing EIRENE runtime for JONAH') 
             !            eirtime = eirtime * 0.25 
          else if (machine2(1:len_trim(machine2)) == 'claire') then 
             CALL WN('GetInput','Increasing EIRENE runtime for CLAIRE') 
             eirtime = eirtime * 2.0 
          else if (machine2(1:len_trim(machine2)) == 'juelich') then 
             CALL WN('GetInput','Increasing EIRENE runtime for JUELICH') 
             eirtime = eirtime * 1.2 
          else if (machine2(1:len_trim(machine2)) == 'joshua') then 
          ELSE 
             CALL WN('GetInput','Unidentified computer: '// &
                  machine2(1:LEN_TRIM(machine2))// &
                  ', EIRENE runtime not modified') 
          end if
       end if
    else if (tag(1:3) == '021') then 
       CALL ReadI(line,eirdata,0,1,'EIRENE input file') 
    else if (tag(1:3) == '022') then 
       CALL ReadI(line,eirmat1,1,5,'EIRENE target material') 
       !      ELSEIF (TAG(1:3).EQ.'023') THEN 
    else if (tag(1:3) == '024') then 
       CALL ReadI(line,eirmat2,1,5,'EIRENE wall material') 
    else if (tag(1:3) == '058') then 
       !...    Load data for EIRENE pressure gauge specifications (additional 
       !       surface data): 
       if (line(7:9) == '1.1') then 
          CALL RdRarn(eirpgdat,eirnpgdat,MAXNAS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,6,'EIRENE pressure gauge specs',ierr) 
          if (ierr /= 0) call er('getinput','eirpgdat',*99) 
          !...      Calculate volume (only toroidal volume at the moment): 
          WRITE(fp,*) 
          WRITE(fp,'(A)') 'Pressure gauge data:' 
          WRITE(fp,'(5X,A4,5A8,2A10)') &
               'No.','x','y','T (deg)','r (m)','Len (m)','P (mTorr)', &
               'Vol (m3)' 
          DO i1 = 1, eirnpgdat 
             !            eirpgdat(i1,8) = PI * eirpgdat(i1,5)**2.0 * 2.0 * PI * 
             !     .                       eirpgdat(i1,2) 
             WRITE(fp,'(5X,I4,5F8.3,1P,2E10.2,0P)') INT(eirpgdat(i1,1)), &
                  (eirpgdat(i1,i2),i2 = 2,8) 
          end do
       ELSE 
          CALL RdRarn(eirpgdat,eirnpgdat,MAXNAS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,4,'EIRENE pressure gauge specs',ierr) 
          if (ierr /= 0) call er('getinput','eirpgdat',*99) 

          DO i1 = 1, eirnpgdat 
             eirpgdat(i1,7) = eirpgdat(i1,5) 
             eirpgdat(i1,5) = eirpgdat(i1,4) 
             eirpgdat(i1,4) = -1.0 
             eirpgdat(i1,6) = -1.0 
             eirpgdat(i1,8) =  1.0 
          end do
          !...      Output gauge specifications: 
          WRITE(fp,*) 
          WRITE(fp,'(A)') 'Pressure gauge data (no volume):' 
          WRITE(fp,'(5X,A4,5A8,2A10)') &
               'No.','x','y','T (deg)','r (m)','Len (m)','P (mTorr)', &
               'Vol (m3)' 
          DO i1 = 1, eirnpgdat 
             WRITE(fp,'(5X,I4,5F8.3,1P,2E10.2,0P)') INT(eirpgdat(i1,1)), &
                  (eirpgdat(i1,i2),i2 = 2,8) 
          end do

       end if

    else if (tag(1:3) == '076') then 
       !...    EIRENE surface properties data (default override): 
       if     (line(7:9) == '1.0') then 

          CALL RdRarn(eirspdat,eirnspdat,MAXNAS3,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,6,'EIRENE surface properties',ierr) 
          if (ierr /= 0) call er('getinput','eirspdat',*99) 

          WRITE(fp,*) 
          WRITE(fp,*) 'SURFACE PROPERTIES:' 
          DO i1 = 1, eirnspdat 
             !...        Assign default value to the surface recycling coefficient: 
             eirspdat(i1,8) = 1.0 

             WRITE(fp,'(A,I4,8F12.6)') '> ',i1,(eirspdat(i1,i2),i2 = 1,8) 
          end do

          if (eirnspdat > 0) eirspdatmode = 1 

       else if (line(7:9) == '2.0'.or.line(7:9) == '2.1') then 

          CALL RdRarn(eirspdat,eirnspdat,MAXNAS3,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,7,'EIRENE surface properties',ierr) 
          if (ierr /= 0) call er('getinput','eirspdat',*99) 

          if (eirnspdat > 0) then 
             if (line(7:9) == '2.0') then 
                eirspdatmode = 2 
             ELSE 
                eirspdatmode = 3 
             end if
          end if

          WRITE(fp,*) 
          WRITE(fp,*) 'SURFACE PROPERTIES:' 
          DO i1 = 1, eirnspdat 
             WRITE(fp,'(A,I4,8F12.6)') '> ',i1,(eirspdat(i1,i2),i2 = 1,8) 
          end do

       else if (line(7:9) == '2.2') then 
          CALL RdRarn(eirspdat,eirnspdat,MAXNAS3,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,8,'EIRENE surface properties',ierr) 
          if (ierr /= 0) call er('getinput','eirspdat',*99) 
          eirspdatmode = 3 
          WRITE(fp,*) 
          WRITE(fp,*) 'SURFACE PROPERTIES:' 
          DO i1 = 1, eirnspdat 
             !...       Move data to accomodate 2.3: 
             eirspdat(i1,10) = eirspdat(i1,9) 
             eirspdat(i1,9) = 1.0 
             WRITE(fp,'(A,I4,9F12.6)') '> ',i1,(eirspdat(i1,i2),i2 = 1,9) 
          end do
       else if (line(7:9) == '2.3') then 
          CALL RdRarn(eirspdat,eirnspdat,MAXNAS3,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,8,'EIRENE surface properties',ierr) 
          if (ierr /= 0) call er('getinput','eirspdat',*99) 
          eirspdatmode = 3 
          WRITE(fp,*) 
          WRITE(fp,*) 'SURFACE PROPERTIES:' 
          DO i1 = 1, eirnspdat 
             WRITE(fp,'(A,I4,9F12.6)') '> ',i1,(eirspdat(i1,i2),i2 = 1,9) 
          end do
       ELSE 
          CALL ER('RUI','Unsupported version for 076 EIRSPDAT',*99) 
       end if

    else if (tag(1:3) == '077') then 
       !...    Additional surfaces for EIRENE 
       if     (line(7:9) == '2.0') then 
          CALL RdRarn(eirasdat,eirnasdat,MAXNAS2,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,9,'EIRENE additional surfaces',ierr) 
          if (ierr /= 0) call er('getinput','eirasdat',*99) 

          eirnasdat = MAX(1,eirnasdat)  ! To preserve the default settings in InitializeUnstructuredInput 
          ! in unstructured_input.f 

          !...      Shift data around so that it is compatible with the existing code.  The 
          !         z-coordinate data is now in elements 8 and 9: 
          DO i1 = 1, eirnasdat 
             if (eirasdat(i1,1) == 1.0.or. &
                  eirasdat(i1,1) == 98.0.or. &
                  eirasdat(i1,1) < 0.0) then 
                z1 = eirasdat(i1,5) 
                DO i2 = 5, 7 
                   eirasdat(i1,i2) = eirasdat(i1,i2+1) 
                end do
                eirasdat(i1,8) = z1 
             end if
          end do
          WRITE(SLOUT,*) 'Additional surface data 2.0:' 
          DO i1 = 1, eirnasdat 
             WRITE(SLOUT,'(I6,10F10.4)') i1,(eirasdat(i1,i2),i2 = 1,10) 
          end do
          ! For deletion... 
       else if (line(7:9) == '1.1') then 
          CALL WN('RUI','Additional wall data format 1.1 soon to be '// &
               'obsolete (unstructured input tag *077)') 
          CALL RdRarn(eirasdat,eirnasdat,MAXNAS2,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,8,'EIRENE additional surfaces',ierr) 
          if (ierr /= 0) call er('getinput','eirasdat',*99) 
          !...      Shift data around so that it is compatible with the existing code.  The 
          !         z-coordinate data is now in elements 8 and 9: 
          DO i1 = 1, eirnasdat 
             z1 = eirasdat(i1,5) 
             DO i2 = 5, 7 
                eirasdat(i1,i2) = eirasdat(i1,i2+1) 
             end do
             eirasdat(i1,8) = z1 
          end do
          WRITE(SLOUT,*) 'Additional surface data:' 
          DO i1 = 1, eirnasdat 
             WRITE(SLOUT,'(I6,9F10.4)') i1,(eirasdat(i1,i2),i2 = 1,9) 
          end do
          !        ELSEIF (line(7:9).EQ.'1.0') THEN 
          !          CALL RdRarn(eirasdat,eirnasdat,MAXNAS2,-MACHHI,MACHHI,.FALSE., 
          !     .              -MACHHI,MACHHI,6,'EIRENE additional surfaces',ierr) 
          !          IF (ierr.NE.0) CALL ER('GetInput','EIRASDAT',*99) 
       ELSE 
          CALL ER('RUI','Unsupported version for 077 EIRASDAT',*99) 
       end if

    else if (tag(1:3) == '012') then 
       CALL ReadI(line,fradd,0,100,'Cells to add to front of ring') 
       ! ---------------------------------------------------------------------- 
    else if (tag(1:3) == '013') then 
       CALL ReadI(line,bkadd,0,100,'Cells to add to end of ring') 
       ! ---------------------------------------------------------------------- 
       ! ---------------------------------------------------------------------- 
    else if (tag(1:3) == '015') then 
       CALL ReadI(line,cmodopt,0,1,'Grid source option') 

       if (cmodopt == 1) thesis = .true. 
       ! ----------------------------------------------------------------------- 
    else if (tag(1:3) == '016') then 
       CALL ReadI(line,stagopt,0,3,'Stagger grid') 
       ! ----------------------------------------------------------------------- 
    else if (tag(1:3) == '017') then 
       CALL ReadI(line,stopopt ,0,999,'Special function setting') 
    else if (tag(1:3) == '063') then 
       CALL ReadI(line,stopopt2,0,999,'Special function setting 2') 
    else if (tag(1:3) == '064') then 
       CALL ReadI(line,stopopt3,0,999,'Special function setting 3') 
    else if (tag(1:3) == '066') then 
       CALL ReadI(line,iflexopt(4),0,999,'Integer flex option 4') 
       if (iflexopt(4) == 30) osm_mode = 2 
    else if (tag(1:3) == '071') then 
       CALL ReadI(line,iflexopt(5),0,999,'Integer flex option 5') 
    else if (tag(1:3) == '073') then 
       CALL ReadI(line,iflexopt(6),0,999,'Integer flex option 6') 
    else if (tag(1:3) == '074') then 
       CALL ReadI(line,iflexopt(7),0,999,'Integer flex option 7') 
    else if (tag(1:3) == '075') then 
       CALL ReadI(line,iflexopt(8),0,999,'Integer flex option 8') 
    else if (tag(1:3) == '065') then 
       CALL ReadR(line,rflexopt(1),-99.0,1.0,'Real flex option 1') 
    else if (tag(1:3) == '067') then 
       CALL ReadR(line,rflexopt(2),-1.0,99.0,'Real flex option 2') 
    else if (tag(1:3) == '068') then 
       CALL ReadR(line,rflexopt(3),0.0,1.0E+8,'Real flex option 3') 
    else if (tag(1:3) == '069') then 
       CALL ReadR(line,rflexopt(4),-1.0,100.0,'Real flex option 4') 
    else if (tag(1:3) == '070') then 
       CALL ReadR(line,rflexopt(5),0.0,10.0,'Real flex option 5') 
    else if (tag(1:3) == '072') then 
       CALL ReadR(line,rflexopt(6),0.0,1.0E+25,'Real flex option 6') 
       ! ----------------------------------------------------------------------- 

       ! ---------------------------------------------------------------------- 
    else if (tag(1:3) == '025') then 
       CALL ReadI(line,haldata,0,1,'CMOD Halpha data') 
       ! ----------------------------------------------------------------------- 
    else if (tag(1:3) == '026') then 
       CALL ReadI(line,pincode,0,5,'PIN neutral code option') 
       ! ----------------------------------------------------------------------- 
    else if (tag(1:3) == '027') then 
       CALL ReadI(line,tarsource,0,8,'Target data source') 
    else if (tag(1:3) == '078') then 
       CALL Read2R(line,tarshift(IKHI),tarshift(IKLO),-HI,HI,'T-shift') 
    else if (tag(1:3) == '085') then 
       CALL ReadI(line,tarshiftopt,0,1,'Target shift option') 
    else if (tag(1:3) == '081') then 
       !...    High index target data for boundary condition relaxation: 
       if     (line(7:9) == '1.1') then 
          CALL RDRARN(LPDATI2,NLPDATI2,MAXINS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,6,'INNER RELAX. DATA',IERR) 
          if (ierr /= 0) call er('readunstructuredinput', &
               '081 error',*99) 
       ELSE 
          CALL ER('ReadUnstructuredInput','Unknown 081 version',*99) 
       end if
    else if (tag(1:3) == '082') then 
       !...    Low index target data for boundary condition relaxation: 
       if     (line(7:9) == '1.1') then 
          CALL RDRARN(LPDATO2,NLPDATO2,MAXINS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,6,'OUTER RELAX. DATA',IERR) 
          if (ierr /= 0) call er('readunstructuredinput', &
               '082 error',*99) 
       ELSE 
          CALL ER('ReadUnstructuredInput','Unknown 082 version',*99) 
       end if

    else if (tag(1:3) == '083') then 
       !...    Low index target data for SOL21,24 over-ride: 
       if     (line(7:9) == '1.1') then 
          CALL RDRARN(s21_datao,s21_ndatao,MAXNRS,-MACHHI,MACHHI, &
               .FALSE.,-MACHHI,MACHHI,9,'OUTER SOL21 DATA',IERR) 
          if (ierr /= 0) call er('readunstructuredinput', &
               '083 error',*99) 
       else if (line(7:9) == '1.0') then 
          CALL RDRARN(s21_datao,s21_ndatao,MAXNRS,-MACHHI,MACHHI, &
               .FALSE.,-MACHHI,MACHHI,7,'OUTER SOL21 DATA',IERR) 
          if (ierr /= 0) call er('readunstructuredinput', &
               '083 error',*99) 
          DO i1 = 1, s21_ndatao 
             s21_datao(i1,9 ) = -1.0 
             s21_datao(i1,10) = -1.0 
          end do
       ELSE 
          CALL ER('ReadUnstructuredInput','Unknown 083 version',*99) 
       end if
    else if (tag(1:3) == '084') then 
       !...    High index target data for SOL21,24 over-ride: 
       if     (line(7:9) == '1.1') then 
          CALL RDRARN(s21_datai,s21_ndatai,MAXNRS,-MACHHI,MACHHI, &
               .FALSE.,-MACHHI,MACHHI,9,'INNER SOL21 DATA',IERR) 
          if (ierr /= 0) call er('readunstructuredinput', &
               '084 error',*99) 
       else if (line(7:9) == '1.0') then 
          CALL RDRARN(s21_datai,s21_ndatai,MAXNRS,-MACHHI,MACHHI, &
               .FALSE.,-MACHHI,MACHHI,7,'INNER SOL21 DATA',IERR) 
          if (ierr /= 0) call er('readunstructuredinput', &
               '084 error',*99) 
          DO i1 = 1, s21_ndatai 
             s21_datai(i1,9 ) = -1.0 
             s21_datai(i1,10) = -1.0 
          end do
       ELSE 
          CALL ER('ReadUnstructuredInput','Unknown 084 version',*99) 
       end if

    else if (tag(1:3) == '088') then 
       !...    High index target data for interpolation: 
       version = GetInputVersion(line) 
       if     (version == 1.0) then 
          tarintermode(IKHI) = 1.0 
          CALL RDRARN(tarinter(1,1,IKHI),tarninter(IKHI), &
               MAXNRS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,6, &
               'HIGH INDEX INTER DATA',IERR) 
          if (ierr /= 0) call er('readunstructuredinput', &
               '088 error',*99) 
          DO i1 = 1, tarninter(IKHI) 
             WRITE(SLOUT,'(A,I6,1P,4E10.2,0P)') &
                  '088:',i1,(tarinter(i1,i2,IKHI),i2 = 1,4) 
          end do
       ELSE 
          tarintermode(IKHI) = 0.0 
          CALL RDRARN(tarinter(1,1,IKHI),tarninter(IKHI), &
               MAXNRS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,3, &
               'HIGH INDEX INTER DATA',IERR) 
          if (ierr /= 0) call er('readunstructuredinput', &
               '088 error',*99) 

          DO i1 = 1, tarninter(IKHI) 
             WRITE(SLOUT,'(A,I6,1P,4E10.2,0P)') &
                  '088:',i1,(tarinter(i1,i2,IKHI),i2 = 1,4) 
          end do
       end if

    else if (tag(1:3) == '089') then 
       !...    Low index target data for interpolation: 
       if     (line(7:9) == '1.0') then 
          tarintermode(IKLO) = 1.0 
          CALL RDRARN(tarinter(1,1,IKLO),tarninter(IKLO), &
               MAXNRS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,6, &
               'LOW INDEX INTER DATA',IERR) 
          if (ierr /= 0) call er('readunstructuredinput', &
               '089 error',*99) 
       ELSE 
          tarintermode(IKLO) = 0.0 
          CALL RDRARN(tarinter(1,1,IKLO),tarninter(IKLO), &
               MAXNRS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,3, &
               'LOW INDEX INTER DATA',IERR) 
          if (ierr /= 0) call er('readunstructuredinput', &
               '089 error',*99) 
       end if
       ! ----------------------------------------------------------------------- 
    else if (tag(1:3) == '028') then 
       CALL ReadR(line,grd_minpl,LO,HI,'Min poloidal side length') 
       ! ----------------------------------------------------------------------- 
    else if (tag(1:3) == '029') then 
       CALL ReadI(line,grd_refine,0,13,'Grid refinement option') 
       ! ----------------------------------------------------------------------- 

       ! ----------------------------------------------------------------------- 
    else if (tag(1:3) == '032') then 
       CALL ReadR(line,grd_range,-HI,HI,'Reginement region') 
       !      ELSEIF (TAG(1:3).EQ.'G01') THEN 
       !...    Not sure if Dave has reserved G01 already, and the web server is down: 
       !        CALL ReadR(line,grd_thresh,0.0,100.0,'Refinement threshold') 
       ! ----------------------------------------------------------------------- 
    else if (tag(1:3) == '033') then 
       CALL ReadR(line,osm_range,0.0,1.0,'SOL 22 error region') 
    else if (tag(1:3) == '042') then 
       CALL ReadR(line,switch(SWPOW2),-1.0,22.0,'SOL 22 pow dist') 
    else if (tag(1:3) == '059') then 
       CALL ReadR(line,switch(SWPOW3),-1.0,19.0,'SOL 22 mock pow dist') 
    else if (tag(1:3) == '043') then 
       CALL ReadR(line,switch(SWION2),-1.0,21.0,'SOL 22 ion dist') 
    else if (tag(1:3) == '044') then 
       CALL ReadI(line,osm_store,-1,100,'load stored PIN sources') 
    else if (tag(1:3) == '046') then 
       CALL ReadI(line,osm_probe,-2,5,'Probe for pressure reference') 
    else if (tag(1:3) == '048') then 
       CALL ReadI(line,s21_mode,0,5,'PIN source use option') 
    else if (tag(1:3) == '050') then 
       CALL ReadI(line,osm_matchs,0,2,'Match plasma at symmetry point') 
    else if (tag(1:3) == '051') then 
       CALL ReadI(line,osm_matchp,0,4,'Match pressure at probe loc.') 
    else if (tag(1:3) == '053') then 
       CALL ReadI(line,osm_symopt,0,5,'SOL 22 symmetry point option') 
    else if (tag(1:3) == '055') then 
       CALL ReadI(line,osm_matcht,0,6,'Match temp. at probe loc.') 
    else if (tag(1:3) == '060') then 
       CALL ReadI(line,osm_preopt,-1,9,'Use SOL22p prescription') 
    else if (tag(1:3) == '061') then 
       CALL ReadI(line,osm_recopt,0,7,'Calculate recombination adj') 
    else if (tag(1:3) == '087') then 
       CALL ReadI(line,osmmock,0,1,'Mock power term option') 
    else if (tag(1:3) == '079') then 
       !...    Plasma data for uniform private flux zone: 
       if (line(7:9) == '1.0') then 
          CALL RdRarn(osmppv,osmnppv,MAXNRS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,5,'Uniform PP data',ierr) 
          if (ierr /= 0) call er('getinput','osmppv',*99) 
          !           DO i1 = 1, osmnppv 
          !            WRITE(0,'(A,1P,6E10.2,0P)') '> ',(osmppv(i1,i2),i2=1,6) 
          !          ENDDO 
       ELSE 
          CALL ER('RUI','Unsupported version for 070 OSMPPV',*99) 
       end if

       ! ----------------------------------------------------------------------- 
    else if (tag(1:3) == '030') then 
       CALL ReadI(line,rel_opt,0,3,'Relaxation option') 
    else if (tag(1:3) == '031') then 
       CALL ReadR(line,rel_frac,0.0,1.0,'PIN relaxation fraction') 
    else if (tag(1:3) == '086') then 
       CALL ReadI(line,relreset,0,2,'Reset PIN sources with each step') 
    else if (tag(1:3) == '034') then 
       CALL ReadI(line,rel_nstep,1,100,'Number of steps') 
    else if (tag(1:3) == '035') then 
       CALL ReadI(line,rel_niter,1,100,'Iteratons per step') 
    else if (tag(1:3) == '036') then 
       CALL ReadR(line,rel_pace,-1.0,9.0,'Target relaxation pace') 
    else if (tag(1:3) == '045') then 
       CALL Read2R(line,rel_bound1,rel_bound2,0.0,1.0,'BC relax. bnd') 
    else if (tag(1:3) == '047') then 
       CALL ReadR(line,rel_tol,0.01,100.0,'Pressure balance tolerance') 
    else if (tag(1:3) == '052') then 
       CALL ReadI(line,osm_powopt,0,2,'Mock power option') 
    else if (tag(1:3) == '041') then 
       CALL Read2I(line,osm_watch1,osm_watch2,0,MAXNRS,'SOL 22 watch') 
    else if (tag(1:3) == '062') then 

       !       Load data for ionisation front specification for SOL22p: 

       CALL RdRarn(osm_ionfnt,osm_nfnt,MAXFNT,-MACHHI,MACHHI, &
            .FALSE.,-MACHHI,MACHHI,2,'Ionisation front spec.     ',ierr) 

       if (ierr /= 0) call er('getinput','reading osm_ionfnt',*99) 
    else if (tag(1:3) == '080') then 
       !...    Relaxation over-ride: 
       if     (line(7:9) == '1.0') then 
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE., &
               0.0,MACHHI,6,'Relaxation data 1.0',IERR) 
          relmode = 0 
          DO i1 = 1, rel_ndata 
             WRITE(fp,'(A,1P,7E10.2)') 'RD: ',(rel_data(i1,i2),i2 = 1,7) 
          end do
       else if (line(7:9) == '1.1') then 
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,8,'Relaxation data 1.1',IERR) 
          relmode = 1 
          DO i1 = 1, rel_ndata 
             WRITE(fp,'(A,1P,9E10.2)') 'RD: ',(rel_data(i1,i2),i2 = 1,9) 
          end do
       else if (line(7:9) == '1.2') then 
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,3,'Relaxation data 1.2',IERR) 
          relmode = 2 
          DO i1 = 1, rel_ndata 
             WRITE(fp,'(A,1P,9E10.2)') 'RD: ',(rel_data(i1,i2),i2 = 1,4) 
          end do
       else if (line(7:9) == '1.3') then 
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,8,'Relaxation data 1.3',IERR) 
          relmode = 3 
          DO i1 = 1, rel_ndata 
             WRITE(fp,'(A,1P,9E10.2)') 'RD: ',(rel_data(i1,i2),i2 = 1,9) 
          end do
       else if (line(7:9) == '1.4') then 
          !...      Over-riding the additional cell plasma parameters: 
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,2,'Relaxation data 1.4',IERR) 
          relmode = 4 
          DO i1 = 1, rel_ndata 
             WRITE(fp,'(A,1P,9E10.2)') 'RD: ',(rel_data(i1,i2),i2 = 1,2) 
          end do
       else if (line(7:9) == '1.5') then 
          !...      Over-riding the additional cell plasma parameters and bulk plasma values: 
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,6,'Relaxation data 1.5',IERR) 
          relmode = 5 
          DO i1 = 1, rel_ndata 
             WRITE(fp,'(A,1P,9E10.2)') 'RD: ',(rel_data(i1,i2),i2 = 1,7) 
          end do
       else if (line(7:9) == '1.6') then 
          !...      Over-riding puff strength: 
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,1,'Relaxation data 1.6',IERR) 
          relmode = 6 
          DO i1 = 1, rel_ndata 
             WRITE(fp,'(A,1P,2E10.2)') 'RD: ',(rel_data(i1,i2),i2 = 1,2) 
          end do
       else if (line(7:9) == '1.7') then 
          !...      Over-riding bulk plasma values: 
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,4,'Relaxation data 1.7',IERR) 
          relmode = 7 
          DO i1 = 1, rel_ndata 
             WRITE(fp,'(A,1P,5E10.2)') 'RD: ',(rel_data(i1,i2),i2 = 1,5) 
          end do
       else if (line(7:9) == '2.0') then 
          !...      Flexible over-ride option: 
          CALL RDRARN(rel_data,rel_ndata,MAXINS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,6,'Relaxation data 2.0',IERR) 
          relmode = 20 
          WRITE(fp,*) 'RELMODE:',relmode 
          DO i1 = 1, rel_ndata 
             WRITE(fp,'(A,1P,7E10.2)') 'RD: ',(rel_data(i1,i2),i2 = 1,7) 
          end do
       ELSE 
          CALL ER('RUI','Unsupported version for 080 REL_DATA',*99) 
       end if
       !...    Set RELMODE to zero if no data is listed in the input file: 
       if (rel_ndata == 0) relmode = 0 
       ! ---------------------------------------------------------------------- 
    else if (tag(1:3) == '037') then 
       CALL ReadI(line,adp_opt,0,2,'Grid adaptation option') 
    else if (tag(1:3) == '038') then 
       CALL ReadI(line,adp_region,1,3,'Adaptation region') 
    else if (tag(1:3) == '039') then 
       CALL ReadR(line,adp_upper,1.0,HI,'Bound for refinement') 
    else if (tag(1:3) == '040') then 
       CALL ReadR(line,adp_lower,1.0,HI,'Bound for coarsification') 
       ! ----------------------------------------------------------------------- 
    else if (tag(1:3) == '056') then 
       CALL ReadI(line,prb_align,0,1,'Adjust position of FSP data') 
    else if (tag(1:3) == '057') then 
       CALL ReadR(line,prb_shift,-5.0,99.0,'Forced FSP shift') 
       if (prb_shift /= 99.0) prb_shift = prb_shift / 1.0e+03 

       ! ----------------------------------------------------------------------- 

       ! ...next...090 -- but it is in use! 


       ! ----------------------------------------------------------------------- 
    else if (tag(1:3) == '090') then 
       CALL ReadI(line,outtarget,0,1,'Output target data') 
    else if (tag(1:3) == '091') then 
       CALL ReadI(line,outwall,0,1,'Output wall data') 
    else if (tag(1:3) == '092') then 
       CALL ReadI(line,outtarget,0,1,'Output target data') 
    else if (tag(1:3) == '093') then 
       CALL ReadI(line,outcell,0,1,'Output cell data') 
    else if (tag(1:3) == '094') then 
       CALL ReadI(line,outgeom,0,1,'Output geometry data') 
    else if (tag(1:3) == '095') then 
       CALL ReadI(line,outpoly,0,1,'Output polygon data') 
    else if (tag(1:3) == '096') then 
       CALL ReadI(line,outpin,0,1,'Output PIN data') 
    else if (tag(1:3) == '097') then 
       CALL ReadI(line,out_source,0,1,'PIN sources') 
    else if (tag(1:3) == '098') then 
       CALL ReadI(line,out_plasma,0,1,'PIN plasma') 
    else if (tag(1:3) == '099') then 
       CALL ReadI(line,out_geom  ,0,1,'PIN geometry') 


       !===================================================
       !
       ! TAG Series E
       !
       !===================================================


    else if (tag(1:3) == 'E11') then 
       CALL ReadI(line,eirbgk    ,0,5,'EIRENE n-n collision option') 
    else if (tag(1:3) == 'E12') then 
       CALL ReadI(line,eiropacity,-5,6,'EIRENE Lyman alpha opacity') 
    else if (tag(1:3) == 'E13') then 
       CALL ReadI(line,eirnstrata,0,0,'EIRENE strata specification') 
    else if (tag(1:3) == 'E14') then 
       CALL ReadI(line,eircxd2   ,0,2,'EIRENE CX D2+ production') 
    else if (tag(1:3) == 'E15') then 
       CALL ReadI(line,eirph2    ,0,1,'EIRENE p-H2 collisions inc.') 
    else if (tag(1:3) == 'E16') then 
       !...    EIRENE puffing surface data: 
       if     (line(7:9) == '1.0') then 
          CALL RdRarn(eirpuff,eirnpuff,MAXNAS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,7,'EIRENE puffing surfaces',ierr) 
          if (ierr /= 0) call er('getinput','E16 eirpuff 1',*99) 
          if (eirnpuff > 0) eirpmode = 1 
          DO i1 = 1, eirnpuff 
             WRITE(fp,'(A,I4,6F10.5)') '> ',i1,(eirpuff(i1,i2),i2 = 1,6) 
          end do
       else if (line(7:9) == '1.1') then 
          CALL RdRarn(eirpuff,eirnpuff,MAXNAS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,5,'EIRENE puffing surfaces',ierr) 
          if (eirnpuff > 0) eirpmode = 2 
          if (ierr /= 0) call er('getinput','E16 eirpuff 2',*99) 
          DO i1 = 1, eirnpuff 
             WRITE(fp,'(A,I4,6F10.5)') '> ',i1,(eirpuff(i1,i2),i2 = 1,6) 
          end do
       else if (line(7:9) == '1.2') then 
          CALL RdRarn(eirpuff,eirnpuff,MAXNAS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,7,'EIRENE puffing surfaces',ierr) 
          if (eirnpuff > 0) eirpmode = 3 
          if (ierr /= 0) call er('getinput','E16 eirpuff 2',*99) 
          DO i1 = 1, eirnpuff 
             WRITE(fp,'(A,I4,8F10.5)') '> ',i1,(eirpuff(i1,i2),i2 = 1,8) 
          end do
       else if (line(7:9) == '1.3') then 
          CALL RdRarn(eirpuff,eirnpuff,MAXNAS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,7,'EIRENE puffing surfaces',ierr) 
          if (eirnpuff > 0) eirpmode = 4 
          if (ierr /= 0) call er('getinput','E16 eirpuff 2',*99) 
          DO i1 = 1, eirnpuff 
             WRITE(fp,'(A,I4,8F10.5)') '> ',i1,(eirpuff(i1,i2),i2 = 1,8) 
          end do
       else if (line(7:9) == '1.4') then 
          CALL RdRarn(eirpuff,eirnpuff,MAXNAS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,9,'EIRENE puffing surfaces',ierr) 
          if (eirnpuff > 0) eirpmode = 4 
          if (ierr /= 0) call er('getinput','E16 eirpuff 4',*99) 
          DO i1 = 1, eirnpuff 
             eirpuff (i1,11) = 0.0 
             WRITE(fp,'(A,I4,10F10.5)') '> ',i1,(eirpuff(i1,i2),i2 = 1,10) 
          end do
       else if (line(7:9) == '1.5') then 
          CALL RdRarn(eirpuff,eirnpuff,MAXNAS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,10,'EIRENE puffing surfaces',ierr) 
          if (eirnpuff > 0) eirpmode = 4 
          if (ierr /= 0) call er('getinput','E16 eirpuff 1.5',*99) 
          DO i1 = 1, eirnpuff 
             WRITE(fp,'(A,I4,11F10.5)') '> ',i1,(eirpuff(i1,i2),i2 = 1,11) 
          end do
       else if (line(7:9) == '2.0') then 
          READ(5,*) cdum1,eirnpuff 
          DO i1 = 1, eirnpuff 
             READ(5,*) eirpuff(1:14,i1),eircpuff(i1) 
             WRITE(fp,'(A,I4,14E10.2,1X,A10)') 'puff: ',i1, &
                  (eirpuff(i2,i1),i2 = 1,14),eircpuff(i1) 
          end do
       ELSE 
          CALL ER('RUI','Unsupported version for E16 EIRPUFF',*99) 
       end if
    else if (tag(1:3) == 'E17') then 
       CALL ReadI(line,eirniter,0,100,'EIRENE self-iterations') 
    else if (tag(1:3) == 'E18') then 
       CALL ReadR(line,eirrinteg,-1.0,1.0,'EIRENE surface par reflec') 
    else if (tag(1:3) == 'E19') then 
       CALL ReadR(line,eireinteg,-1.0,1.0,'EIRENE surface eng reflec') 
    else if (tag(1:3) == 'E20') then 
       CALL ReadR(line,eirermin,-1.0,10.0,'EIRENE thermal cut off') 
    else if (tag(1:3) == 'E21') then 
       CALL Read2I(line,eirtrc1,eirtrc2,0,999999,'Particle trck range') 
    else if (tag(1:3) == 'E22') then 
       CALL ReadR(line,eirzaa,-10.0,1000.0,'toroidal circumference') 
    else if (tag(1:3) == 'E23') then 
       if     (line(7:9) == '1.1') then 
          CALL RdIarn(eiraout,eirnaout,MAXASD,1,2,.TRUE., &
               -99,MAXASCDAT,5,'EIRENE add cell output',ierr) 
          if (ierr /= 0) call er('getinput','E23 eiraout 1',*99) 

          if (eirnaout == 0) then 
             WRITE(fp,*) 'SETTING DEFAULT FOR EIRAOUT' 
             eirnaout = 1 
             eiraout(1,1) =  2 
             eiraout(1,2) =  0 
             eiraout(1,3) = -MAXASCDAT 
             eiraout(1,4) =  MAXASCDAT 
             eiraout(1,5) =  0 
             eiraout(1,6) =  0 
          end if

          WRITE(fp,*) 'ADDITIONAL CELL OUTPUT LIST 1.1:' 
          DO i1 = 1, eirnaout 
             WRITE(fp,'(A,I4,6I10)') '> ',i1,(eiraout(i1,i2),i2 = 1,6) 
          end do
       else if (line(7:9) == '1.0') then 
          CALL RdIarn(eiraout,eirnaout,MAXASD,1,2,.TRUE., &
               -99,MAXASCDAT,4,'EIRENE add cell output',ierr) 
          if (ierr /= 0) call er('getinput','E23 eiraout 1',*99) 

          !...      Shift data to account for version 1.1 format: 
          DO i1 = 1, eirnaout 
             DO i2 = 6,3,-1 
                eiraout(i1,i2) = eiraout(i1,i2-1) 
             end do
             eiraout(i1,2) = 1 
          end do

          if (eirnaout == 0) then 
             WRITE(fp,*) 'SETTING DEFAULT FOR EIRAOUT' 
             eirnaout = 1 
             eiraout(1,1) =  2 
             eiraout(1,2) = -MAXASCDAT 
             eiraout(1,3) =  MAXASCDAT 
             eiraout(1,4) =  0 
             eiraout(1,5) =  0 
          end if

          WRITE(fp,*) 'ADDITIONAL CELL OUTPUT LIST 1.0:' 
          DO i1 = 1, eirnaout 
             WRITE(fp,'(A,I4,6I10)') '> ',i1,(eiraout(i1,i2),i2 = 1,6) 
          end do
       ELSE 
          CALL ER('ReadUnstructuredInput','Invalid version',*99) 
       end if
    else if (tag(1:3) == 'E24') then 
       CALL ReadR(line,eirtorfrac,0.0001,1.0,'toroidal fraction') 
    else if (tag(1:3) == 'E25') then 
       CALL ReadR(line,eirsrcmul,0.0,1.0E+6,'source multiplication') 
    else if (tag(1:3) == 'E26') then 
       CALL ReadI(line,eirfuji,-1,1,'Fujimoto D2+ rates') 
    else if (tag(1:3) == 'E27') then 
       CALL ReadR(line,eiralloc,0.0,1.0,'CPU-time weighting') 
    else if (tag(1:3) == 'E28') then 
       if     (line(7:9) == '1.0') then 
          CALL RdRarn(eirstrata,eirnstrata,MAXSTRATA,-MACHHI,MACHHI, &
               .FALSE.,-MACHHI,MACHHI,4,'EIRENE stratum data',ierr) 
          if (ierr /= 0) call er('getinput','E28 eirstrata 1',*99) 
          WRITE(fp,*) 'EIRENE STRATUM DATA 1.0:' 
          DO i1 = 1, eirnstrata 
             WRITE(fp,'(A,I4,5F10.5)') '> ',i1,(eirstrata(i1,i2),i2 = 1,5) 
          end do
       ELSE 
          CALL ER('ReadUnstructuredInput','Invalid E28 version',*99) 
       end if
    else if (tag(1:3) == 'E29') then 
       CALL ReadI(line,eiracx    ,0,1,'EIRENE CX reaction included') 
    else if (tag(1:3) == 'E30') then 
       CALL ReadI(line,eird2ion  ,0,1,'EIRENE D2 ionisation') 
    else if (tag(1:3) == 'E31') then 
       CALL ReadI(line,eird2dis  ,0,1,'EIRENE D2 dissociation') 
    else if (tag(1:3) == 'E32') then 
       if     (line(7:9) == '1.0') then 
          call clearbuf
          CALL RdRarn(eiriontime,eirniontime,MAXIONTIME,-MACHHI,MACHHI, &
               .FALSE.,-MACHHI,MACHHI,7,'EIRENE time-to-ion',ierr) 
          if (ierr /= 0) call er('getinput','E32 eiriontime',*99) 
          WRITE(fp,*) 'EIRENE TIME-TO-IONISATION PARAMETERS 1.0:' 
          DO i1 = 1, eirniontime 
             WRITE(fp,'(A,I4,8F10.5)') '> ',i1,(eiriontime(i1,i2),i2 = 1,8) 
             if (eiriontime(i1,5) > maxbin-2) &
                  CALL ER('ReadUnstructuredInput','Number of time bins '// &
                  'requested in E32 exceeds maximum',*99) 
          end do
       ELSE 
          CALL ER('ReadUnstructuredInput','Invalid E32 version',*99) 
       end if
    else if (tag(1:3) == 'E33') then 
       !...    Some custom input here: 
       if     (line(7:9) == '1.0') then 
          READ(5,*) cdum1,eirntally 
          if (eirntally > maxtally) &
               CALL ER('ReadUnstructuredInput','Too many volume tallies '// &
               'specified.  Increase MAXTALLY.',*99) 
          DO i1 = 1, eirntally 
             READ(5,*) (eirtally(i1,i2),i2 = 1,5) 
             WRITE(PINOUT,'(I2,4(1X,A))') i1, &
                  (eirtally(i1,i2)(1:LEN_TRIM(eirtally(i1,i2))),i2 = 1,5) 
          end do
       ELSE 
          CALL ER('ReadUnstructuredInput','Invalid E32 version',*99) 
       end if
    else if (tag(1:3) == 'E34') then 
       !...    Load EIRENE atomic data multipliers: 
       !         1 - D ionisation 
       !         2 - D CX 
       !         3 - D-D collisions 
       !         4 - D-D2 collisions 
       !         5 - D2 ionisation 
       !         6 - D2 dissociation 1 
       !         7 - D2 dissociation 2 
       !         8 - D2-D collisions 
       !         9 - D2-D2 collisions 
       !        10 - charge exchange production of D2 
       !        11 - D+ recombination rate 
       !        12 - p-D2 collision rate 

       CALL ReadIR(line,idum1,rdum1,1,12,'EIRENE atomic data mul.') 
       eirscale(idum1) = rdum1 
    else if (tag(1:3) == 'E35') then 
       CALL ReadI(line,eirnsection,1,100,'Number of toroidal sections') 
    else if (tag(1:3) == 'E36') then 
       !...    Specify regions where EIRENE non-standard default surfaces are 
       !       transparent: 
       if     (line(7:9) == '1.0') then 
          CALL RdRarn(eirtrans,eirntrans,MAXTOR,-MACHHI,MACHHI, &
               .FALSE.,-MACHHI,MACHHI,2,'EIRENE time-to-ion',ierr) 
          if (ierr /= 0) call er('getinput','E36 eirtrans',*99) 
          WRITE(fp,*) 
          WRITE(fp,*) 'EIRENE TRANSPARENT SURFACE SPAN 1.0:' 
          DO i1 = 1, eirntrans 
             WRITE(fp,'(A,I4,3F10.5)') '> ',i1,(eirtrans(i1,i2),i2 = 1,3) 
          end do
       ELSE 
          CALL ER('ReadUnstructuredInput','Invalid E36 version',*99) 
       end if
    else if (tag(1:3) == 'E37') then 
       CALL ReadI(line,eirntorseg,0,500,'Num sections in EIRENE appro') 
    else if (tag(1:3) == 'E38') then 
       CALL ReadR(line,eirdtimv,0.0,100.0,'Time dependent mode int.') 
    else if (tag(1:3) == 'E39') then 
       CALL ReadI(line,eirtrim,0,1,'EIRENE TRIM database option') 
    else if (tag(1:3) == 'E40') then 
       IF (.TRUE.) THEN 
          CALL RDRARN(eirtri,eirntri,MAXNRS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,4,'Additional triangles',IERR) 
          DO i1 = 1, eirntri 
             WRITE(fp,'(A,1P,1P,8E10.2,0P)') 'TRI:',(eirtri(i1,i2),i2 = 1,8) 
          end do
       end if

    else if (tag(1:3) == 'E41') then 
       CALL ReadI(line,eirphoton,-1,2,'photons being followed') 
       ! ---------------------------------------------------------------------- 







       
       !===================================================
       !
       ! TAG Series A
       !
       !===================================================
       !

    elseif (tag(1:3) .eq. 'A01') then
       !  Inputs for TAG series A
       !
       !   Sample input 
       !   '+A01 Ref DIV Title' 'Case name and header'
       call divrd(TITLE, 'TITLE FOR RUN', IERR)
    elseif (tag(1:3) .eq. 'A02') then
       !   Sample input 
       !   '+A02 DIV Desc:' 'Text Describing the case'
       call divrd(desc, 'DESCRIPTION OF RUN', IERR)
    elseif (tag(1:3) .eq. 'A03') then
       !   Sample input 
       !   '+A03 Equil File Name' 'sonnet_oskn_c08'
       call divrd(equil,'equilibrium file name',ierr)
    elseif (tag(1:3) .eq. 'A04') then
       !   Sample input 
       !   '+A04    Print option  (0 reduced, 1 full)               '      0
       call divrd(CPRINT,.true., 0 ,.true.,  19 ,'PRINT OPTION',IERR)

    else if (tag(1:3) == 'A05') then 

       !     TAG A05: 

       !       ne_opt - selects how electron density is defined 
       !       0 -> ne = nb    1 -> ne = nb + sigma nz (from fluid code) 

       CALL ReadI(line,ne_opt,0,1,'Electron density option') 

       ! ----------------------------------------------------------------------- 


       !     TAG A06: 

    else if (tag(1:3) == 'A06') then 

       !       Option to write a JET TRAN file for POST-PROCESSOR use 
       !       Written to Unit 41 - 0 = off - 1 = 0n .. default ON 


       CALL ReadI(line,write_tran,0,1,'TRAN FILE PRINT OPTION') 

       ! ----------------------------------------------------------------------- 


       !     TAG A07: 

    else if (tag(1:3) == 'A07') then 

       !       Option to write a netcdf version of the raw data file 
       !       0 = off   1 = 0n .. default OFF 

       CALL ReadI(line,netcdf_opt,0,1,'WRITE NETCDF FORMAT RAW OUTPUT') 



       !===================================================
       !
       ! TAG Series B
       !
       !===================================================
       !
       !  Inputs for TAG series B
       !
    elseif (tag(1:3) .eq. 'B01') then
       !   Sample input 
       !   '+B01    Special plasma parameter           Rspec        '     1
       call divrd(IRSPEC,.FALSE., 0 ,.FALSE., 0 ,'SPEC PLASMA RSPEC',IERR)
    elseif (tag(1:3) .eq. 'B02') then
       !   Sample input 
       !   '+B02    For average density "near" target  xnear (m)    '    0.6
       call divrd(CXNEAR,.TRUE. ,0.0,.FALSE.,0.0,'NEAR PARAM  XNEAR ',IERR)
    elseif (tag(1:3) .eq. 'B03') then
       !   Sample input 
       !   '+B03    For average density "near" target  ynear (m) +/-'    0.6
       call divrd(CYNEAR,.TRUE. ,0.0,.FALSE.,0.0,'NEAR PARAM  YNEAR ',IERR)
    elseif (tag(1:3) .eq. 'B04') then
       !   Sample input 
       !   '+B04    Debug atoms (0 off, >0 print every nth timestep)'      0
       call divrd(ISTEPN ,.FALSE., 0 ,.FALSE., 0 ,'*** DEBUG NEUT ***',IERR)
    elseif (tag(1:3) .eq. 'B05') then
       !   Sample input 
       !   '+B05    Debug ions  (0 off, >0 print every nth timestep)'      0
       call divrd(ISTEPL ,.FALSE., 0 ,.FALSE., 0 ,'*** DEBUG DIV  ***',IERR)
    elseif (tag(1:3) .eq. 'B06') then
       !   Sample input 
       !   '+B06    Debug ion velocity  (0 off, >1 on)              '      0
       call divrd(ISTEPV ,.FALSE., 0 ,.FALSE., 0 ,'**DEBUG VELOCITY**',IERR)
    elseif (tag(1:3) .eq. 'B07') then
       !   Sample input 
       !   '+B07 TN442 Z-value defining divertor region  (m)         '   1.7
       call divrd(CZD   ,.FALSE.,0.0,.FALSE.,0.0, 'Z DIVERTOR LIMIT',IERR)
    elseif (tag(1:3) .eq. 'B08') then
       !   Sample input 
       !   '+B08 TN    Ring number for detailed background data      '    0
       call divrd(CIRHR,  .TRUE., 0,.TRUE.,MAXNRS,'HI-RES RING     ',IERR)




       !===================================================
       !
       ! TAG Series C
       !
       !===================================================
       !
       !  Inputs for TAG series C
       !
    elseif (tag(1:3) .eq. 'C01') then
       !   Sample input 
       !   '+C01 ' 'Set of S-distances for ion leakage diagnostic(m)'
       !   'TN982     Number of S-values  :-'  0
       ! Example
       !'TN982     Number of S-values  :-'  4
       !   5.0
       !   10.0
       !   15.0
       !   20.0
       call divrd(cleaks,cleaksn,maxpts,0.0,machhi,.TRUE.,'LEAK- S',IERR)
    elseif (tag(1:3) .eq. 'C02') then
       !   Sample input 
       !   '+C02 TN1303 Dperp Extractor - methods used - 0 1 2       '     2
       call divrd(dpmethod,.TRUE.,0  ,.true.,2 ,'DPERP EXT METHOD',IERR)
    elseif (tag(1:3) .eq. 'C03') then
       !   Sample input 
       !   '+C03 TN1310 Dperp Ext. 0=only to Xpoint 1=Full Field Line'     1
       call divrd(dpsuml,.TRUE.,0  ,.true.,1 ,'DPERP EXT SUM LIMIT',IERR)
    elseif (tag(1:3) .eq. 'C04') then
       !   Sample input 
       !   '+C04 TN1310 Dperp Ext. Outer Ring Losses    0=off 1=on   '     1
       call divrd(dpouter,.TRUE.,0  ,.true.,1 ,'DP EXT OUTER RING',IERR)
    elseif (tag(1:3) .eq. 'C05') then
       !   Sample input 
       !   '+C05 TN1309 Dperp Ext. Dperp  Convection    0=off 1=on   '     1
       call divrd(dpconv,.TRUE.,0  ,.true.,1 ,'DP EXT CONVECT LOSS',IERR)
    elseif (tag(1:3) .eq. 'C06') then
       !   Sample input 
       !   '+C06 TN1311 Dperp Ext. 1/2 Cell Flux Corr.  0=off 1=on   '     0
       call divrd(dpfluxopt,.TRUE.,0  ,.true.,1 ,'CELL CENTRE FLUX',IERR)
    elseif (tag(1:3) .eq. 'C07') then
       !   Sample input 
       !   '+C07 TN1311 Dperp Ext. Calc. Average Dperp. 0=off N=on   '     1
       call divrd(dpavopt,.TRUE.,0  ,.false.,1 ,'AVERAGE DPERP OPT',IERR)
    elseif (tag(1:3) .eq. 'C08') then
       !   Sample input 
       !   '+C08 TN1311 Dperp Ext. Major Radius Corr.   0=off 1=on   '     0
       call divrd(dprcopt,.TRUE.,0  ,.true.,2 ,'MAJOR RADIUS CORR',IERR)
    elseif (tag(1:3) .eq. 'C09') then
       !   Sample input 
       !   '+C09 TN1314 Dperp Ext. Gradient Smoothing   0=off N=on   '     0
       call divrd(dpsmooth,.true.,0,.false.,0 ,'GRADIENT SMOOTHING',IERR)
    elseif (tag(1:3) .eq. 'C10') then
       !   Sample input 
       !   '+C10 TN     Dperp Ext. Gradient Calc Meth     -1,0,1     '     0
       call divrd(dpnav,.true.,-1,.true.,2 ,   'GRADIENT CALC METH',IERR)
    elseif (tag(1:3) .eq. 'C11') then
       !   Sample input 
       !   '+C11 TN     Dperp Ext. Cross-field Area  0-centre 1-bound'     0
       call divrd(dparea,.TRUE.,0  ,.true.,1 ,'CELL BOUND AREA',IERR)
    elseif (tag(1:3) .eq. 'C12') then
       !   Sample input 
       !   '+C12 TN     Dperp Ext. Power Loss Terms     0-off 1-on   '     1
       call divrd(dpploss,.TRUE.,0  ,.true.,2 ,'POWER LOSS TERMS',IERR)
    elseif (tag(1:3) .eq. 'C13') then
       !   Sample input 
       !   '+C13 TN     Dperp Ext. Non-ortho Correction 0-off 1-on   '     1
       call divrd(dporth,.TRUE.,0 ,.true.,1,'DP NON-ORTH GRADIENT',IERR)
    elseif (tag(1:3) .eq. 'C14') then
       !   Sample input 
       !   '+C14 TN     Dperp Ext. Pei Correction Factor             '    1.0
       call divrd(dppei,.TRUE.,0.0,.false.,0.0,'PEI Correction',IERR)
    elseif (tag(1:3) .eq. 'C15') then
       !   Sample input 
       !   '+C15 TN1373 Dperp Ext. Source Recycle Correction Factor  '    1.0
       call divrd(dprec,.TRUE.,0.0,.false.,0.0,'EXT Recycle Frac',ierr)
    elseif (tag(1:3) .eq. 'C16') then
       !   Sample input 
       !   '+C16 TN1445 Dperp Ext. Dperp/Xperp Fixed Ratio  0.0=off  '    0.0
       call divrd(dpxpratio,.FALSE.,0.0,.false.,0.0,'Dp/Xp Ratio',ierr)
    elseif (tag(1:3) .eq. 'C17') then
       !   Sample input 
       !   '+C17 Vertical   Reciprocating Probe - Intersection #  '     1
       !
       !     Reciprocating/Fast Scanning Probe - R location.
       !
       !     .... and horizontal probe location
       !
       call divrd(rlocnum,.true.,1,.false.,0, 'Crossing number',ierr)
    elseif (tag(1:3) .eq. 'C18') then
       !   Sample input 
       !   '+C18 Vertical   Reciprocating Probe - R-Value         '    1.94
       call divrd(crploc,.FALSE.,0.0,.false.,0.0,'Rec.probe loc',ierr)
    elseif (tag(1:3) .eq. 'C19') then
       !   Sample input 
       !   '+C19   Horizontal Reciprocating Probe - Intersection #  '     1
       call divrd(zlocnum,.true.,1,.false.,0, 'Crossing number',ierr)
    elseif (tag(1:3) .eq. 'C20') then
       !   Sample input 
       !   '+C20   Horizontal Reciprocating Probe - Z-Value         '    0.0
       call divrd(czploc,.FALSE.,0.0,.false.,0.0,'Rec.probe loc',ierr)

    else if (tag(1:3) == 'C21') then 
       ! ----------------------------------------------------------------------- 

       !     TAG C21 - read in a myultiplier for PINQE as used in the 
       !                  Dperp/Xperp extractor 

       CALL ReadR(line,dp_pinqe_mult,0.0,HI,'PINQE Multiplier for Extractor') 

    else if (tag(1:3) == 'C22') then 
       ! ----------------------------------------------------------------------- 

       !     TAG C22 - Calculate a velocity shifted line profile for the 
       !               specified neutral impurity line. 
       CALL ReadI(line,line_profile_opt,0,1,'Line Profile Calculation Switch') 

       !       Read in the additional required data. 

       !       Read in ADAS Selector data 
       call rdg1(line3,lp_ADASID,lp_ADASYR,lp_ADASEX,lp_ISELE,lp_ISELR,lp_ISELX,lp_ISELD,IERR) 

       !        write(0,'(a,4i5)') 'ADAS:',lp_isele, 
       !     >          len(lp_adasid),len(lp_adasex), 
       !     >          lp_adasyr 


       !        write(0,'(a,4i5)') 'ADAS:',lp_isele, 
       !     >          len(lp_adasid),len(lp_adasex), 
       !     >          lp_adasyr 
       !        write(0,'(a,a,a))') 'ADAS:',lp_adasid,':' 
       !        write(0,'(a,a,a))') 'ADAS:',lp_adasex,':' 


       !       Read in LOS data. 

       call rd_lp_los(line3,lp_robs,lp_zobs,lp_theta,lp_dtheta,lp_instrument_width,lp_bin_width,ierr) 

       !       Convert angles in degrees to radians 

       lp_theta = lp_theta * degrad 
       lp_dtheta = lp_dtheta * degrad 




       !===================================================
       !
       ! TAG Series D
       !
       !===================================================
       !
       !  Inputs for TAG series D
       !
    elseif (tag(1:3) .eq. 'D01') then
       !   Sample input 
       !   '+D01    Source Data Option 0-Nocorona 1-Adas            '     1
       call divrd(cdatopt,.true.,0,.true.,3, 'Rad/ioniz data source',ierr)
    elseif (tag(1:3) .eq. 'D02') then
       !   Sample input 
       !   '+D02    USERID for ADAS H data (*=use central database) '   '*'
       call divrd(useridh,'ADAS H userid',ierr)
    elseif (tag(1:3) .eq. 'D03') then
       !   Sample input 
       !   '+D03    Year for ADAS H data                            '    96
       call divrd(iyearh,.true., 0,.true.,99,'ADAS H year          ',ierr)
    elseif (tag(1:3) .eq. 'D04') then
       !   Sample input 
       !   '+D04    USERID for ADAS Z data (*=use central database) '   '*'
       call divrd(useridz,'ADAS Z userid',ierr)
    elseif (tag(1:3) .eq. 'D05') then
       !   Sample input 
       !   '+D05    Year for ADAS Z data                            '    96
       call divrd(iyearz,.true., 0,.true.,99,'ADAS Z year          ',ierr)
    elseif (tag(1:3) .eq. 'D06') then
       !   Sample input 
       !   '+D06 B2FRATES name:' '/home/david/divimp/adpak/C_rates.strahl'
       call divrd(mcfile,'Name of File for ADPAK/INEL atomic data',ierr)
    elseif (tag(1:3) .eq. 'D07') then
       !   Sample input 
       !   '+D07    Sputter data option 1-old 2-93                  '     6
       call divrd(CSPUTOPT,.TRUE., 1,.TRUE., 8,'SPUTTER SOURCE OPT ',IERR)
    elseif (tag(1:3) .eq. 'D08') then
       !   Sample input 
       !   '+D08    Chemical Sputter Data Option                    '    11
       call divrd(CCHEMOPT,.TRUE., 1,.TRUE.,12,'CHEMSPUT SOURCE OPT',IERR)
    elseif (tag(1:3) .eq. 'D09') then
       !   Sample input 
       !   '+D09 TN1481 BG Ion Mom.Tran.Coll.Coef.     (kelighi)     '    5.0E-16
       call divrd(kelighi,.true.,0.0,.false.,0.0,'MTC COEFF 1 Imp-I',ierr)
    elseif (tag(1:3) .eq. 'D10') then
       !   Sample input 
       !   '+D10 TN1481 BG Neutral Mom.Tran.Coll.Coef. (kelighg)     '    5.0E-16
       call divrd(kelighg,.true.,0.0,.false.,0.0,'MTC COEFF 2 Imp-N',ierr)
    elseif (tag(1:3) .eq. 'D11') then
       !   Sample input 
       !   '+D11    Characteristic energy Ebd          Ebd   (eV)   '    0.0
       call divrd(CEBD,  .TRUE. ,0.0,.FALSE.,0.0,'CHAR ENERGY EBD',   IERR)
    elseif (tag(1:3) .eq. 'D12') then
       !   Sample input 
       !   '+D12    Neutral hydrogen density parameter Nhc   (m**-3)'    1.0E15
       call divrd(CNHC,  .FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  NHC   ',IERR)
    elseif (tag(1:3) .eq. 'D13') then
       !   Sample input 
       !   '+D13                                       Nho   (m**-3)'    3.0E18
       call divrd(CNHO,  .FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  NHO   ',IERR)
    elseif (tag(1:3) .eq. 'D14') then
       !   Sample input 
       !   '+D14                                       lamhx (m)    '    0.02
       call divrd(CLAMHX,.FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  LAMHX ',IERR)
    elseif (tag(1:3) .eq. 'D15') then
       !   Sample input 
       !   '+D15                                       lamhy (m)    '    0.11
       call divrd(CLAMHY,.FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  LAMHY ',IERR)
    elseif (tag(1:3) .eq. 'D16') then
       !   Sample input 
       !   '+D16    Constant for CX Recomb option 2    Vcx   (m/s)  '    2.4E4
       call divrd(CVCX  ,.TRUE., 0.0,.FALSE.,0.0,'CXREC CONSTANT VCX',IERR)
    elseif (tag(1:3) .eq. 'D17') then
       !   Sample input 
       !   '+D17    Threshold yield for self-sputtering      (eV)   '    0.1
       call divrd(CTRESH,.FALSE.,0.0,.FALSE.,0.0,'SELF-SPU THRESHOLD',IERR)
    elseif (tag(1:3) .eq. 'D18') then
       !   Sample input 
       !   '+D18    Bombarding ion charge state        Zimp         '      2
       call divrd(CBOMBZ,.TRUE.,  0 ,.FALSE., 0 ,'BOMBION CHARGE    ',IERR)
    elseif (tag(1:3) .eq. 'D19') then
       !   Sample input 
       !   '+D19    Bombion type 0Zb 1H 2D 3T 4He4 5C 6Zi 7O        '      5
       call divrd(CBOMBF,.TRUE.,  0 ,.TRUE.,  7 ,'BOMBION FLAG 0:7  ',IERR)
    elseif (tag(1:3) .eq. 'D20') then
       !   Sample input 
       !   '+D20    Ionisation rate factor for neutrals          IRF'    1.0
       call divrd(CIRF  ,.TRUE., 0.0,.FALSE.,0.0,'IONISE RATE FACTOR',IERR)
    elseif (tag(1:3) .eq. 'D21') then
       !   Sample input 
       !   '+D21    Sputtering Enhancement Factor                SEF'    1.0
       call divrd(CSEF  ,.TRUE., 0.0,.FALSE.,0.0,'SPUT ENHANCE FACT.',IERR)
    elseif (tag(1:3) .eq. 'D22') then
       !   Sample input 
       !   '+D22 ' 'Set of Yield Modifiers for Primary, Secondary neutrals'
       !   '      Number of rows of (X,Mpt,Mst,Mct,Mpw,Mcw,Refl) data-'  0
       !         0.0   0.0    1.0     1.0    1.0    1.0    1.0     1.0
       !
       !---- READ IN YIELD MODIFIER FUNCTIONS
       !
       ! note default values are set so that only modifications need to be read in
       ! if you want to change the defaults applied to all elements then the first line in cymfs should be for
       ! indices 0.0  0.0 
       call divrd(CYMFS,NYMFS,MAXPTS+1,-machhi,machhi,.true.,-machhi,machhi,7,'SET YIELD(P,S,CT,PW,CW) VALS',IERR)
    elseif (tag(1:3) .eq. 'D23') then
       !   Sample input 
       !   '+D23 TN? Fixed Yield Value for Sputter Data Option 4     '     0.001
       call divrd(const_yield,.true.,0.0,.FALSE.,0.0,'Fixed Yield',IERR)
    elseif (tag(1:3) .eq. 'D24') then
       !   Sample input 
       !   '+D24 TN1209 Target Temperature (K) for Chem. Sputt. Opt. '   300.0
       call divrd(ctargt ,.true.,0.0,.FALSE.,0.0,'Target Temp (K) ',IERR)
    elseif (tag(1:3) .eq. 'D25') then
       !   Sample input 
       !   '+D25       Main Wall Temperature (K) for Chem. Sputt.   '   300.0
       call divrd(cwallt ,.true.,0.0,.FALSE.,0.0,'Wall Temp (K)   ',IERR)
    elseif (tag(1:3) .eq. 'D26') then
       !   Sample input 
       !   '+D26 TN1450 PP Wall   Temperature (K) for Chem. Sputt.   '   300.0
       call divrd(cwallp ,.true.,0.0,.FALSE.,0.0,'PP Wall Temp (K)',IERR)
    elseif (tag(1:3) .eq. 'D27') then
       !   Sample input 
       !   '+D27 ' 'TN1450 Wall Temperatures (K) for specific segments'
       !   '      Number of segment ranges (Index1 Index2 Temp):'      0
       call divrd(walltemp,NWLTEMP,MAXPTS,-MACHHI,MACHHI,.FALSE.,0.0,MACHHI,2,'WALL SEGMENT TEMPERATURES',IERR)
    elseif (tag(1:3) .eq. 'D28') then
       !   Sample input 
       !   '+D28    Temperature Gradient Coefficient parameter  ZO  '      4
       call divrd(CZO   ,.TRUE. , 0 ,.FALSE., 0 ,'ZO PARAMETER      ',IERR)
    elseif (tag(1:3) .eq. 'D29') then
       !   Sample input 
       !   '+D29      CEMAXF factor for sputtering velocities       '   1.0
       call divrd(CEMAXF,.FALSE.,0.0,.FALSE.,0.0,'EMAX-FACTOR',      IERR)
    elseif (tag(1:3) .eq. 'D30') then
       !   Sample input 
       !   '+D30 TN521 Impact Energy for wall launch Vel. dist (eV)  '  100.0
       call divrd(CEIMP, .TRUE., 0.0,.FALSE.,0.0,'EIMP- WALL LAUNCH',IERR)
    elseif (tag(1:3) .eq. 'D31') then
       !   Sample input 
       !   '+D31 TN83? Maximum Number of sputtered generations       '   75
       call divrd(CMAXGENS,.TRUE., 0,.FALSE.,  0, 'MAX. GENERATIONS',IERR)
    elseif (tag(1:3) .eq. 'D32') then
       !   Sample input 
       !   '+D32 TN1007 ABSFAC or Power - Specified - Use if > 0.0   '   0.0
       call divrd(nabsfac,.false.,0.0,.false.,0.0, 'ABSFAC modifier',ierr)
    elseif (tag(1:3) .eq. 'D33') then
       !   Sample input 
       !   '+D33  TN1200 Stgrad - Distance where Tgrad forces -> 0    '   1.0
       call divrd(Cstgrad,.TRUE.,0.0,.TRUE.,1.0, 'TGRAD FORCES -> 0',IERR)
    elseif (tag(1:3) .eq. 'D34') then
       !   Sample input 
       !   '+D34    H Recombination Calculation Option 4-Adas 0+-oth'     4
       call divrd(crecopt,.TRUE.,0  ,.true. ,4  ,'Recomb. Calc Opt',IERR)
    elseif (tag(1:3) .eq. 'D35') then
       !   Sample input 
       !   '+D35    Recombination Limit Cut-OFF Temperature (eV)    '    0.0
       call divrd(treccut,.TRUE.,0.0,.false.,0.0,'Recomb. cutoff T',IERR)
    elseif (tag(1:3) .eq. 'D36') then
       !   Sample input 
       !   '+D36 TN1429 T-GRAD Force Modification - Factor Applied   '    1.0
       !
       !     Factor for temperature gradient force modifier option
       !
       call divrd(fgradfact,.false.,0.0,.false.,0.0,'FGRAD MOD FACT',ierr)

    else if (tag(1:3) == 'D37') then 
       ! ----------------------------------------------------------------------- 

       !     TAG D37 and D38 
       !     ADAS IONIZATION AND RECOMBINATION RATE MODIFIERS 

       !     These values can be used to modify the rates read in from ADAS 
       !     The default values should always be set to 1.0 and these 
       !     should be used with care if used at all. 

       !     TAG D37 - adas_iz_rate_mult 
       CALL ReadR(line,adas_iz_rate_mult,0.0,HI, &
            'ADAS Ionization rate multiplier') 

    else if (tag(1:3) == 'D38') then 
       !     TAG D38 - adas_rec_rate_mult 

       CALL ReadR(line,adas_rec_rate_mult,0.0,HI, &
            'ADAS Recombination rate multiplier') 


    else if (tag(1:3) == 'D39') then 
       ! ----------------------------------------------------------------------- 

       !     TAG D39 : Alternate Sputter data specifier - used to select one of 
       !               several custom sputter datasets - usually based 
       !               on different impact angles 

       !     Secondary sputter data specifier 

       CALL ReadR(line,extra_sputter_angle,-10.0,90.0,'Extra Sputter Angle Opt') 

    else if (tag(1:3) == 'D40') then 
       ! ----------------------------------------------------------------------- 

       !     TAG D40 : Flux fraction for alternate bombarding ion sputter 
       !               calculations. Goes with CBOMBF and CBOMBZ to allow for 
       !               specification of trace impurity sputtering in cases 
       !               where hydrogenic sputtering is expected to be negligible. 

       !               Ideally this should be expanded to allow for 
       !               hydrogen+a specifiable distribution of impurity charge states 


       !     Bombarding ion flux fraction 

       CALL ReadR(line,cbomb_frac,0.0,1.0,'Bombarding ion flux fraction') 

    else if (tag(1:3) == 'D41') then 
       !     TAG D41-43 
       !     These options are related to sputtering from SiC using the SiC 
       !     mixed-material model in Abrams NF 2022. Note this requires 
       !     csputopt = 8. 
       !     D41: Usage switch. Useful for comparisons of SiC to its 
       !          constituents using the same yield data that is hardocded in 
       !          for SiC. 
       !          1 = Normal. Target is SiC. 
       !          2 = Target is graphite. 
       !          3 = Target is silicon. 
       !     D42: Flux fraction of the carbon in the incoming plasma flux. 
       !     D43: Flux fraction of the silicon in the incoming plasma flux. 

       call readi(line, mm_usage, 0, 2, &
            'SiC mixed-material model usage switch') 
    else if (tag(1:3) == 'D42') then 
       call readr(line, frac_c, 0.0, 1.0, &
            'Fraction of C in SiC mixed-material model') 
    else if (tag(1:3) == 'D43') then 
       call readr(line, frac_si, 0.0, 1.0, &
            'Fraction of Si in SiC mixed-material model') 

    else if (tag(1:3) == 'D44') then 
       !     TAG D44 
       !     Usage switch for TiB2 (0) or ZrB2 (1). 
       call readi(line, tib2_or_zrb2, 0, 1, &
            'TiB2 or ZrB2 usage switch') 




       !===================================================
       !
       ! TAG Series F
       !
       !===================================================
       !
       !  Inputs for TAG series F
       !
    elseif (tag(1:3) .eq. 'F01') then
       !   Sample input 
       !   '+F01    Read EDGE2D BG for reference  0=No  1=Yes       '     0
       call divrd(cre2d ,.TRUE., 0,.true.,2 ,'READ E2D BG FOR REF  ',ierr)
    elseif (tag(1:3) .eq. 'F02') then
       !   Sample input 
       !   '+F02    Use EDGE2D Target Data Option 0=Off 1=Reg 2=Flux'     0
       call divrd(e2dtargopt,.TRUE., 0,.true.,6 ,'EDGE2D TARG COND',ierr)
    elseif (tag(1:3) .eq. 'F03') then
       !   Sample input 
       !   '+F03    Plasma condition for missing SOL rings     CNIWA'      1
       call divrd(CNIWA ,.TRUE. , 0 ,.FALSE., 1 ,'CNIWA OPTION      ',IERR)
    elseif (tag(1:3) .eq. 'F04') then
       !   Sample input 
       !   '+F04 ' 'EDGE1D/2D Deuterium drift vel. mult. factor VMF '
       !   '            Number of VMF blocks                          '      0
       !   '            Ring Range :-' -20     -30
       !   '            J0 & J1    :-'   5       5
       !   '            VMF0,1,2   :-'   1.000   1.000   1.000
       call divrd('VEL. MULT. FACTOR',IERR)
    elseif (tag(1:3) .eq. 'F05') then
       !   Sample input 
       !   '+F05 SONNET-Number of Fluids in B2 Background Plasma File'     7
       call divrd(nfla,.TRUE.,1,.true.,maxnfla,'NUM.FLUIDS IN B2 BG',IERR)
    elseif (tag(1:3) .eq. 'F06') then
       !   Sample input 
       !   '+F06       Read Aux. Background Input File    0=off 1=on'     0
       call divrd(readaux,.TRUE.,0,.true.,3,'READ AUX BG FILE',IERR)

    else if (tag(1:3) == 'F11') then 
       ! ----------------------------------------------------------------------- 

       !     TAG F11 - Uedge background option 

       CALL ReadI(line,uedge_bg,0,1,'UEDGE background option') 

    else if (tag(1:3) == 'F12') then 

       ! ----------------------------------------------------------------------- 

       !     Options affecting the calculation/interpretation of fluid code 
       !     target conditions. 

       !     Option F12 sets groups of these flags so that input can be 
       !     simplified. Any individual flag options appearing later in the 
       !     input file will overwrite the values specified by the overall 
       !     flag option. 

       CALL ReadI(line,fc_target_calc_option,0,2,'Fluid Code Target Data Calculation Option') 

       !       Assign values to related sub-options 

       !       Base JET standard 

       if (fc_target_calc_option == 0) then 
          fc_v_calc_opt  = 0 
          fc_te_calc_opt = 1 
          fc_ti_calc_opt = 2 
          fc_ne_calc_opt = 2 

          !       Base UEDGE standard 

       else if (fc_target_calc_option == 1) then 
          fc_v_calc_opt  = 0 
          fc_te_calc_opt = 0 
          fc_ti_calc_opt = 0 
          fc_ne_calc_opt = 2 

          !       Base Old Divimp/B2 standard 

       else if (fc_target_calc_option == 2) then 
          fc_v_calc_opt  = 1 
          fc_te_calc_opt = 2 
          fc_ti_calc_opt = 2 
          fc_ne_calc_opt = 2 

          !       Alternate B2/B2.5/Eirene 

       else if (fc_target_calc_option == 3) then 
          fc_v_calc_opt  = 1 
          fc_te_calc_opt = 1 
          fc_ti_calc_opt = 1 
          fc_ne_calc_opt = 2 
       end if

    else if (tag(1:3) == 'F13') then 
       CALL ReadI(line,fc_ne_calc_opt,0,2,'Fluid Code Target Ne Data Calculation Option') 
    else if (tag(1:3) == 'F14') then 
       CALL ReadI(line,fc_te_calc_opt,0,2,'Fluid Code Target Te Data Calculation Option') 
    else if (tag(1:3) == 'F15') then 
       CALL ReadI(line,fc_ti_calc_opt,0,2,'Fluid Code Target Ti Data Calculation Option') 
    else if (tag(1:3) == 'F16') then 
       CALL ReadI(line,fc_v_calc_opt,0,1,'Fluid Code Target Vb Data Calculation Option') 
    else if (tag(1:3) == 'F17') then 
       CALL ReadI(line,fc_v_interp_opt,0,1,'Fluid Code Cell Edge Value Interpretation Option') 
    else if (tag(1:3) == 'F18') then 
       CALL ReadI(line,e2dformopt,0,2,'Fluid Code Plasma File Format Specifier') 

    else if (tag(1:3) == 'F19') then 
       !     For readaux = 3 (read auxiliary fluid code data option 3 - SOLPS4.3) 
       !     fort.44 file containing neutral densities. 

       CALL ReadI(line,e2dneut_select,1,3,'Specify which block of impurity neutral data to read') 


    else if (tag(1:3) == 'F20') then 
       !     F20: 

       !     e2dion_select = 1 

       !     The e2dnzs data stored in fort.31 can contain multiple fluid species 
       !     in addition to H+. This quantity specifies an offset into this data 
       !     so that the code can start reading the correct impurity into the e2d 
       !     fluid code arrays like e2dnzs. This is required in DIVIMP so it can 
       !     include meaningful comparisons between the DIVIMP and fluid code results. 

       CALL ReadI(line,e2dion_select,1,100, &
            'Specify where to start reading the fluid code ion data') 



       !===================================================
       !
       ! TAG Series N
       !
       !===================================================

       !
       !  Inputs for TAG series N
       !
    elseif (tag(1:3) .eq. 'N01') then
       !   Sample input 
       !   '+N01    Launch option  0distrib 1point 2asymp 3tip 4wall'     0
       call divrd(CNEUTB,.TRUE., 0,.TRUE., 7,'LAUNCH OPTION        ',IERR)
    elseif (tag(1:3) .eq. 'N02') then
       !   Sample input 
       !   '+N02    Vel/angle flag 0-11                             '     0
       call divrd(CNEUTC,.TRUE., 0,.TRUE.,20,'VEL/ANGLE FLAG       ',IERR)
    elseif (tag(1:3) .eq. 'N03') then
       !   Sample input 
       !   '+N03 TN487 Supplemental Launch Option (as above)         '     0
       !
       !     THESE VARIABLES (CNEUTH AND CNEUTI) WERE ADDED TO ALLOW LAUNCHES FROM THE WALLS,
       !     IN ADDITION TO THE NORMAL PLATE LAUNCHES. HOWEVER, THE
       !     CAPABILITY CAN BE USED TO SIMULATE TWO DIFFERING DISTRIBUTIONS
       !     OF PARTICLES COMING FROM THE PLATES. PERHAPS DUE TO
       !     PHYSICAL AND CHEMICAL SPUTTERING SIMULATNEOUSLY.
       !
       call divrd(CNEUTH,.TRUE.,-1,.TRUE., 5,'SUP LAUNCH OPTION    ',IERR)
    elseif (tag(1:3) .eq. 'N04') then
       !   Sample input 
       !   '+N04 TN487 Supplemental V/A flag      (as above)         '     3
       call divrd(CNEUTI,.TRUE.,-1,.TRUE.,19,'SUP VEL/ANGLE FLAG   ',IERR)
    elseif (tag(1:3) .eq. 'N05') then
       !   Sample input 
       !   '+N05    Initial Neutral Vel/Ang flag (-1=above,0-13)    '    -1
       call divrd(NVAOPT,.TRUE.,-1,.TRUE.,19,'NEUT VEL/ANGLE FLAG  ',IERR)
    elseif (tag(1:3) .eq. 'N06') then
       !   Sample input 
       !   '+N06 TN1490 Supplemental 2D Neutral Launch 0=off 1=UEDGE '     0
       call divrd(neut2d_opt,.TRUE.,0,.TRUE., 2,'EXTRA 2D LAUNCH OPTION',IERR)
    elseif (tag(1:3) .eq. 'N07') then
       !   Sample input 
       !   '+N07 TN1490 V/A Flag for 2D Neutral Launch (as regular)  '    3
       call divrd(neut2d_vaopt,.TRUE.,-1,.TRUE.,20,'EXTRA 2D V/A FLAG',IERR)
    elseif (tag(1:3) .eq. 'N08') then
       !   Sample input 
       !   '+N08    Sputter option 0std 1special 2mix 3self 4selfva1'     0
       call divrd(CNEUTD,.TRUE., 0,.TRUE., 8,'SPUTTER OPTION       ',IERR)
    elseif (tag(1:3) .eq. 'N09') then
       !   Sample input 
       !   '+N09 TN1209 Secondary Sputter Option                     '     0
       call divrd(CNEUTD2,.TRUE.,-1,.TRUE.,8,'2ND SPUTTER OPTION   ',IERR)
    elseif (tag(1:3) .eq. 'N10') then
       !   Sample input 
       !   '+N10    Normal option  0std 1fromT                      '     0
       call divrd(CNEUTE,.TRUE., 0,.TRUE., 2,'NORMAL OPTION        ',IERR)
    elseif (tag(1:3) .eq. 'N11') then
       !   Sample input 
       !   '+N11    NEUT spreading 0off 1on                         '     0
       call divrd(CNEUTF,.TRUE., 0,.TRUE., 1,'NEUT SPREADING       ',IERR)
    elseif (tag(1:3) .eq. 'N12') then
       !   Sample input 
       !   '+N12 TN1488 Neutral Impurity Velocity Type Option        '     0
       call divrd(cneutvel,.true.,0,.true.,2,'IMP.NEUT.VEL.TYPE.OPT',ierr)
    elseif (tag(1:3) .eq. 'N13') then
       !   Sample input 
       !   '+N13 T   Neutral Wall Reflection 0-off 1-on              '     0
       ! ammod - added options 3 and 4 to NRFOPT
       call divrd(NRFOPT,.TRUE., 0,.TRUE., 4,'NEUTRAL REFLECTION   ',IERR)
    elseif (tag(1:3) .eq. 'N14') then
       !   Sample input 
       !   '+N14 TN1481 Imp Neutral Momentum.Tran.Coll Opt 0=off 1=on'     2
       call divrd(mtcopt,.true.,0,.true.,2,'NEUT.MOM.TRAN.COLL OPT ',ierr)
    elseif (tag(1:3) .eq. 'N15') then
       !   Sample input 
       !   '+N15    Measure theta from T (degrees to X=0) for launch'   90.0
       call divrd(CSNORM,.FALSE.,0.0,.FALSE.,0.0,'MEASURE THETA FROM',IERR)
    elseif (tag(1:3) .eq. 'N16') then
       !   Sample input 
       !   '+N16 ' 'TN487 Launch probability modifiers for each     '
       !   '    TN487 Wall segment range  #1  #2  mod.  :- '           0
       call divrd(WLPROB,NWLPROB,MAXPTS,-MACHHI,MACHHI,.FALSE.,0.0,MACHHI,2,'WALL LAUNCH PROB. MODIFIERS',IERR)
    elseif (tag(1:3) .eq. 'N17') then
       !   Sample input 
       !   '+N17 TN721 Use Wall Probabilities as Absolute 0=No 1=Yes '    0
       call divrd(WLPABS,.TRUE. ,0  ,.TRUE. , 3 ,'ABS WALL PROB    ',IERR)
    elseif (tag(1:3) .eq. 'N18') then
       !   Sample input 
       !   '+N18 A18   Power of cosine release dist. (V/A 12,13)     '   1.0
       call divrd(CNIN,.TRUE.,0.01,.FALSE.,0.0,'COSINE DIST. POWER',IERR)
    elseif (tag(1:3) .eq. 'N19') then
       !   Sample input 
       !   '+N19 TN???? Extra Velocity Multiplier for V/A 14         '   1.0
       call divrd(CVAMULT,.TRUE.,LO,.FALSE.,0.0,'VELOCITY MULTIPLIER',IERR)
    elseif (tag(1:3) .eq. 'N20') then
       !   Sample input 
       !   '+N20 TN???? Velocity Multiplier for Recombined Neutrals  '   1.0
       call divrd(CVRMULT,.TRUE.,LO,.FALSE.,0.0,'REC. VEL. MULT',IERR)

    else if(tag(1:3) == 'N21') then 
       ! ----------------------------------------------------------------------- 

       !     TAG N21 : External sputtering flux data source 
       !               0 = geier file format for Ar 
       !               1 = import divimp charge resolved flux and energy 
       !                   data from a previous divimp run 
       CALL ReadI(line,ext_flx_data_src,0,1, &
            'External sputter flux data source option') 


       !===================================================
       !
       ! TAG Series P
       !
       !===================================================
       !
       !  Inputs for TAG series P
       !
    elseif (tag(1:3) .eq. 'P01') then
       !   Sample input 
       !   '+P01    SOL option     0,1,1a,2,3,4,5,9,10  99file      '    22
       call divrd(CIOPTF,.TRUE.,-1,.TRUE.,99,'SOL OPT              ',IERR)
    elseif (tag(1:3) .eq. 'P02') then
       !   Sample input 
       !   '+P02    Core Option    0,1,2,3                          '     0
       call divrd(ccoreopt,.true.,-1,.true.,28,'Core Option        ',ierr)
    elseif (tag(1:3) .eq. 'P03') then
       !   Sample input 
       !   '+P03    Plasma decay   0std                 99file      '     0
       call divrd(CIOPTG,.TRUE., 0,.TRUE.,99,'PLASMA DECAY OPT     ',IERR)
    elseif (tag(1:3) .eq. 'P04') then
       !   Sample input 
       !   '+P04 ' 'BG PLASMA Options by Ring (PlasDec opts 90&91)  '
       !   '    R1,R2,Sect, PlasDec, Sol, Teg, Tig, Core, Efield'     0
       !
       !     Example
       !     10 23  3.0      4.0 22.0  0.0  0.0   0.0     3.0
       !     26 33  3.0      4.0 22.0  0.0  0.0   0.0     3.0
       call divrd(BGPLASOPT,NBGPLAS,2*MAXNRS,-MACHHI,MACHHI,.FALSE.,0.0,MACHHI,11,'SET OF BG PLASMA OPTIONS BY RING',IERR)
    elseif (tag(1:3) .eq. 'P05') then
       !   Sample input 
       !   '+P05    Trap Tgrad Opt 0off 1on                         '     1
       call divrd(CIOPTO,.TRUE., 0,.TRUE., 4,'TRAP TGRAD OPTION    ',IERR)
    elseif (tag(1:3) .eq. 'P06') then
       !   Sample input 
       !   '+P06    SOL Enhancement Factor - Electric Field    SOLEF'    1.0
       call divrd(CSOLEF,.TRUE., 0.0,.FALSE.,0.0,'SOL  ENHANCE E(S) ',IERR)
    elseif (tag(1:3) .eq. 'P07') then
       !   Sample input 
       !   '+P07    SOL Enhancement Factor - Drift Velocity    SOLVF'    1.0
       call divrd(CSOLVF,.TRUE., 0.0,.FALSE.,0.0,'SOL  ENHANCE VH(S)',IERR)
    elseif (tag(1:3) .eq. 'P08') then
       !   Sample input 
       !   '+P08    SOL1a Factor                               fl   '    0.01
       call divrd(CFL   ,.FALSE.,0.0,.FALSE.,0.0,'SOL1A FACTOR "FL" ',IERR)
    elseif (tag(1:3) .eq. 'P09') then
       !   Sample input 
       !   '+P09    SOL1a Factor                               fs   '    1.0
       call divrd(CFS   ,.FALSE.,0.0,.FALSE.,0.0,'SOL1A FACTOR "FS" ',IERR)
    elseif (tag(1:3) .eq. 'P10') then
       !   Sample input 
       !   '+P10    SOL10 Reversal Mach Number                 fRM  '    1.0
       call divrd(CFRM  ,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACTOR "FRM"',IERR)
    elseif (tag(1:3) .eq. 'P11') then
       !   Sample input 
       !   '+P11    SOL10 factor                               kin  '    1.0
       call divrd(CKIN  ,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACTOR "KIN"',IERR)
    elseif (tag(1:3) .eq. 'P12') then
       !   Sample input 
       !   '+P12    SOL10 factor                               kout '    1.2
       call divrd(CKOUT ,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACTOR"KOUT"',IERR)
    elseif (tag(1:3) .eq. 'P13') then
       !   Sample input 
       !   '+P13    SOL10 factor                               fRmin'    0.01
       call divrd(CFRMIN,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACT "FRMIN"',IERR)
    elseif (tag(1:3) .eq. 'P14') then
       !   Sample input 
       !   '+P14    SOL10 factor                               fRmax'    0.4
       call divrd(CFRMAX,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACT "FRMAX"',IERR)
    elseif (tag(1:3) .eq. 'P15') then
       !   Sample input 
       !   '+P15    SOL6&7 Vb Length factor 1 (* SMAX)         VbL1 '    0.166
       call divrd(CVBL1 ,.TRUE.,0.0 ,.TRUE. ,1.0,'SOL11 LENGTH 1',IERR)
    elseif (tag(1:3) .eq. 'P16') then
       !   Sample input 
       !   '+P16    SOL6&7 Vb multiplication factor 1          VbM1 '    0.0
       call divrd(CVBM1 ,.false.,0.0,.FALSE.,1.0,'SOL11 MULT 1  ',IERR)
    elseif (tag(1:3) .eq. 'P17') then
       !   Sample input 
       !   '+P17    SOL6&7 Vb Length factor 2 (* SMAX)         VbL2 '    0.5
       call divrd(CVBL2 ,.TRUE.,0.0 ,.TRUE. ,1.0,'SOL11 LENGTH 2',IERR)
    elseif (tag(1:3) .eq. 'P18') then
       !   Sample input 
       !   '+P18    SOL6&7 Vb multiplication factor 2          VbM2 '    0.0
       call divrd(CVBM2 ,.false.,0.0,.FALSE.,1.0,'SOL11 MULT 2  ',IERR)
    elseif (tag(1:3) .eq. 'P19') then
       !   Sample input 
       !   '+P19    Power density                      P/A    (W/m2)'    3.0E+07
       call divrd(CPA  ,.FALSE.,0.0,.FALSE.,0.0,'POWER DENSITY', IERR)
    elseif (tag(1:3) .eq. 'P20') then
       !   Sample input 
       !   '+P20    Parallel heat conduction coeff     K0           '    2.0E+03
       call divrd(CK0  ,.FALSE.,0.0,.FALSE.,0.0,'PAR HEAT COND', IERR)
    elseif (tag(1:3) .eq. 'P21') then
       !   Sample input 
       !   '+P21    Parallel heat conduction -ions     K0I          '     58.9
       call divrd(CK0I ,.FALSE.,0.0,.FALSE.,0.0,'PAR HEAT COND IONS',IERR)
    elseif (tag(1:3) .eq. 'P22') then
       !   Sample input 
       !   '+P22    Override input E-field from file E=0  0-off 1-on'      3
       call divrd(OFIELD,.TRUE. , 0 ,.TRUE. , 5 ,'EFIELD=0 FOR FILE',IERR)
    elseif (tag(1:3) .eq. 'P23') then
       !   Sample input 
       !   '+P23 T1433 E-field Opt 4 - Length of E-field region *SMAX'    0.25
       call divrd(CEFLEN,.TRUE. ,0.0,.FALSE.,0.0,'EFIELD SRC LENGTH',IERR)
    elseif (tag(1:3) .eq. 'P24') then
       !   Sample input 
       !   '+P24 T1433 E-field Opt 4 - Te collisionality multiplier  '    2.0
       call divrd(CEFfact,.TRUE. ,0.0,.FALSE.,0.0,'EFIELD COLL-MULT',IERR)
    elseif (tag(1:3) .eq. 'P25') then
       !   Sample input 
       !   '+P25 TN401 Decay length for ionization source Ls  SOL12  '    0.08
       call divrd(CSOLLS,.TRUE.,0.0,.FALSE.,0.0, 'LS DECAY SOL12   ',IERR)
    elseif (tag(1:3) .eq. 'P26') then
       !   Sample input 
       !   '+P26 T     Second decay characteristic length            '    0.08
       call divrd(CSOLLT,.TRUE.,0.0,.FALSE.,0.0, 'SECOND DECAY DIST',IERR)
    elseif (tag(1:3) .eq. 'P27') then
       !   Sample input 
       !   '+P27 T     Source fraction                               '    0.5
       call divrd(CFIZ  ,.TRUE.,0.0,.TRUE. ,1.0, 'SOURCE FRACTION  ',IERR)
    elseif (tag(1:3) .eq. 'P28') then
       !   Sample input 
       !   '+P28 TN401 Decay length for radiative losses  Lr  SOL12  '    0.08
       call divrd(CSOLLR,.TRUE.,0.0,.FALSE.,0.0, 'LR DECAY SOL12   ',IERR)
    elseif (tag(1:3) .eq. 'P29') then
       !   Sample input 
       !   '+P29 TN401 Coefficient for radiative losses   Pr  SOL12  '    1.0
       call divrd(CSOLPR,.TRUE.,0.0,.FALSE.,0.0, 'PR CONST SOL12   ',IERR)
    elseif (tag(1:3) .eq. 'P30') then
       !   Sample input 
       !   '+P30 TN775 Radiation source strength fraction Frr        '    1.0
       call divrd(CSOLFR,.TRUE.,0.0,.FALSE.,0.0, 'FR SOURCE FRAC.  ',IERR)
    elseif (tag(1:3) .eq. 'P31') then
       !   Sample input 
       !   '+P31 TN401 Source Ionization Option 0-lin 1-exp   SOL12  '    1
       call divrd(CSOPT ,.TRUE., 0,.TRUE.,5, 'IONIZATION SOURCE OPT',IERR)
    elseif (tag(1:3) .eq. 'P32') then
       !   Sample input 
       !   '+P32 TN401 Source Radiation Option  0-lin 1-exp   SOL12  '    3
       call divrd(CPOPT ,.TRUE., 0,.TRUE.,3, 'RADIATION SOURCE OPT ',IERR)
    elseif (tag(1:3) .eq. 'P33') then
       !   Sample input 
       !   '+P33 TN660 Imaginary Root Option                  SOL12+ '    1
       call divrd(SROOTOPT,.TRUE.,0,.TRUE.,1,'TREAT NEGATIVE ROOTS ',IERR)
    elseif (tag(1:3) .eq. 'P34') then
       !   Sample input 
       !   '+P34 TN407 Flux Recirculation Option 0-off 1-on          '    0
       call divrd(FLUXROPT,.TRUE.,0,.TRUE.,1,'SRC FLUX RECIRC OPT',  IERR)
    elseif (tag(1:3) .eq. 'P35') then
       !   Sample input 
       !   '+P35 ' 'TN????+407 Set of flux specifications  '
       !   '    Ring , data1(I), data2(O), data3 : Rows - '          0
       call divrd(FLUXINFO,FLUXPTS,MAXINS,-MACHHI,MACHHI,.FALSE.,-machhi,MACHHI,3,'SET OF FLUX DATA',IERR)
    elseif (tag(1:3) .eq. 'P36') then
       !   Sample input 
       !   '+P36 TN408 Calculate SOL iteratively? 0-NO 1-YES         '    0
       call divrd(CITERSOL,.TRUE.,0 ,.TRUE. ,2  ,'DO SOL AFTER PIN ',IERR)
    elseif (tag(1:3) .eq. 'P37') then
       !   Sample input 
       !   '+P37 TN408 Secondary SOL option for iterative calculation'   -2
       !   IF (CSECSOL.EQ.-2) CSECSOL = CIOPTF
       call divrd(CSECSOL ,.TRUE.,-2,.TRUE. ,14 ,'SECONDARY SOL OPT',IERR)
    elseif (tag(1:3) .eq. 'P38') then
       !   Sample input 
       !   '+P38 TN408 Ionization option for iterative SOL           '    3
       call divrd(CSOPT2,.TRUE.,0,.TRUE.,5,  'SECOND ION SOURCE OPT',IERR)
    elseif (tag(1:3) .eq. 'P39') then
       !   Sample input 
       !   '+P39      Number of Pin/SOL iterations                  '    3
       call divrd(NITERSOL,.TRUE.,0 ,.FALSE.,maxpiniter,'NO. OF PIN ITER. ',IERR)
    elseif (tag(1:3) .eq. 'P40') then
       !   Sample input 
       !   '+P40 TN     Te S1 - First  S -value = ctes1 * SMAX       '    0.05
       !
       !     The following values apply to the specifiable SOL option that
       !     allows two part linearly interpolated fits for each of
       !     Te, Ti, ne and vb to be specified.
       !
       !
       ! The following lines specify the parameters for the linearly
       ! interpolated specified BG SOL option. This is purely empirical.
       !                    S-value      Function Value
       ! The form is:         0.0            F0
       !                       S1            F1
       !                       S2            F2
       ! For S>S2 F=F2
       !
       call divrd(ctes1,.TRUE.,0.0,.false.,0.0,'Te distance  1',ierr)
    elseif (tag(1:3) .eq. 'P41') then
       !   Sample input 
       !   '+P41 TN     Te F1 - First  Te-value = ctef1 * te0        '    1.5
       call divrd(ctef1,.TRUE.,0.0,.false.,0.0,'Te fact/mult 1',ierr)
    elseif (tag(1:3) .eq. 'P42') then
       !   Sample input 
       !   '+P42 TN     Te S2 - Second S -value = ctes2 * SMAX       '    0.3
       call divrd(ctes2,.TRUE.,0.0,.false.,0.0,'Te distance  2',ierr)
    elseif (tag(1:3) .eq. 'P43') then
       !   Sample input 
       !   '+P43 TN     Te F2 - Second Te-value = ctef2 * te0        '    2.0
       call divrd(ctef2,.TRUE.,0.0,.false.,0.0,'Te fact/mult 2',ierr)
    elseif (tag(1:3) .eq. 'P44') then
       !   Sample input 
       !   '+P44 TN     Ti S1 - First  S -value = ctis1 * SMAX       '    0.05
       call divrd(ctis1,.TRUE.,0.0,.false.,0.0,'Ti distance  1',ierr)
    elseif (tag(1:3) .eq. 'P45') then
       !   Sample input 
       !   '+P45 TN     Ti F1 - First  Te-value = ctif1 * ti0        '    1.5
       call divrd(ctif1,.TRUE.,0.0,.false.,0.0,'Ti fact/mult 1',ierr)
    elseif (tag(1:3) .eq. 'P46') then
       !   Sample input 
       !   '+P46 TN     Ti S2 - Second S -value = ctis2 * SMAX       '    0.3
       call divrd(ctis2,.TRUE.,0.0,.false.,0.0,'Ti distance  2',ierr)
    elseif (tag(1:3) .eq. 'P47') then
       !   Sample input 
       !   '+P47 TN     Ti F2 - Second Te-value = ctif2 * ti0        '    2.0
       call divrd(ctif2,.TRUE.,0.0,.false.,0.0,'Ti fact/mult 2',ierr)
    elseif (tag(1:3) .eq. 'P48') then
       !   Sample input 
       !   '+P48 TN     Nb S1 - First  S -value = cnes1 * SMAX       '    0.05
       call divrd(cnes1,.TRUE.,0.0,.false.,0.0,'Nb distance  1',ierr)
    elseif (tag(1:3) .eq. 'P49') then
       !   Sample input 
       !   '+P49 TN     Nb F1 - First  Te-value = cnef1 * ne0        '    2.0
       call divrd(cnef1,.TRUE.,0.0,.false.,0.0,'Nb fact/mult 1',ierr)
    elseif (tag(1:3) .eq. 'P50') then
       !   Sample input 
       !   '+P50 TN     Nb S2 - Second S -value = cnes2 * SMAX       '    0.35
       call divrd(cnes2,.TRUE.,0.0,.false.,0.0,'Nb distance  2',ierr)
    elseif (tag(1:3) .eq. 'P51') then
       !   Sample input 
       !   '+P51 TN     Nb F2 - Second Te-value = cnef2 * ne0        '    2.0
       call divrd(cnef2,.TRUE.,0.0,.false.,0.0,'Nb fact/mult 2',ierr)
    elseif (tag(1:3) .eq. 'P52') then
       !   Sample input 
       !   '+P52 TN     vb S1 - First  S -value = cvbs1 * SMAX       '    0.1
       call divrd(cvbs1,.TRUE.,0.0,.false.,0.0,'vb distance  1',ierr)
    elseif (tag(1:3) .eq. 'P53') then
       !   Sample input 
       !   '+P53 TN     vb F1 - First  Te-value = cvbf1 * ti0        '    0.25
       call divrd(cvbf1,.TRUE.,0.0,.false.,0.0,'vb fact/mult 1',ierr)
    elseif (tag(1:3) .eq. 'P54') then
       !   Sample input 
       !   '+P54 TN     vb S2 - Second S -value = cvbs2 * SMAX       '    0.4
       call divrd(cvbs2,.TRUE.,0.0,.false.,0.0,'vb distance  2',ierr)
    elseif (tag(1:3) .eq. 'P55') then
       !   Sample input 
       !   '+P55 TN     vb F2 - Second Te-value = cvbf2 * ti0        '    0.0
       call divrd(cvbf2,.TRUE.,0.0,.false.,0.0,'vb fact/mult 2',ierr)
    elseif (tag(1:3) .eq. 'P56') then
       !   Sample input 
       !   '+P56 TN1424 Core Option4,5- Velocity decay fraction *SMAX'    0.05
       !
       !     Marfe Core options - for option 4,5,6
       !
       call divrd(corefv,.true.,0.0,.true.,0.5,'CORE-VEL FRAC',ierr)
    elseif (tag(1:3) .eq. 'P57') then
       !   Sample input 
       !   '+P57 TN1424 Core Option4,5- Te,Ti    decay fraction *SMAX'    0.05
       call divrd(coreft,.true.,0.0,.true.,0.5,'CORE-TEMP FRAC',ierr)
    elseif (tag(1:3) .eq. 'P58') then
       !   Sample input 
       !   '+P58 TN1427 Core Option4,5- Velocity decay frac 2   *SMAX'    0.5
       call divrd(corefv2,.true.,0.0,.true.,0.5,'CORE-VEL FRAC2',ierr)
    elseif (tag(1:3) .eq. 'P59') then
       !   Sample input 
       !   '+P59 TN1427 Core Option4,5- Te,Ti    decay frac 2   *SMAX'    0.1
       call divrd(coreft2,.true.,0.0,.true.,0.5,'CORE-TEMP FRAC2',ierr)

    else if (tag(1:3) == 'P60') then 
       ! ----------------------------------------------------------------------- 

       !     TAG P60 : Density Gradient Option 
       CALL ReadI(line,ngradopt,0,1,'Density Gradient option') 

    else if (tag(1:3) == 'P61') then 
       ! ----------------------------------------------------------------------- 

       !     TAG P61 : Background Flow Velocity Override Option 

       CALL ReadI(line,override_bg_velocity_opt,0,14, &
            'Background Velocity Override Option') 

    else if (tag(1:3) == 'P62') then 
       ! ----------------------------------------------------------------------- 

       !     TAG P62 : Over-ride Efield target E-field calculation 

       CALL ReadI(line,ofield_targ,1,3, &
            'Override E-field target') 


    else if (tag(1:3) == 'P63') then 
       ! ----------------------------------------------------------------------- 

       !     TAG P63 : External plasma overlay option 
       !               0 = off 
       !               1 = on 
       CALL ReadI(line,external_plasma_overlay,0,1,'External plasma overlay') 


    else if (tag(1:3) == 'P64') then 
       ! ----------------------------------------------------------------------- 

       !     TAG P64 : External plasma overlay file name 
       !               - specifies the name of the file to be loaded 
       !               - full path required unless rundiv script is modified 
       CALL ReadC(line,external_plasma_file,'EXTERNAL PLASMA OVERLAY FILE NAME') 

    else if (tag(1:3) == 'P65') then 
       ! ----------------------------------------------------------------------- 

       !     TAG P65 : SOL option 12/13 etc - additional pressure option 
       !               PMULT - Adds additional pressure to PMULT * PINF 
       !               Over a distance of PDIST * SMAX 

       CALL ReadR(line,sol13_padd,0.0,HI,'SOL13+ ADDITIONAL PRESSURE') 



    else if (tag(1:3) == 'P66') then 
       !----------------------------------------------------------------------- 
       !     TAG P66 : SOL option 12/13 etc - additional pressure option 
       !               PDIST - Adds additional pressure to PMULT * PINF 
       !               Over a distance of PDIST * SMAX 
       CALL ReadR(line,sol13_pdist,0.0,HI,'SOL13+ ADDITIONAL PRESSURE DISTANCE') 



       !===================================================
       !
       ! TAG Series Q
       !
       !===================================================

       !
       !  Inputs for TAG series Q
       !
    elseif (tag(1:3) .eq. 'Q01') then
       !   Sample input 
       !   '+Q01    TeB Gradient   0lin 1lin/lin 2p/a   99file      '     0
       call divrd(CIOPTK,.TRUE., 0,.TRUE.,99,'TEB GRADIENT OPTION  ',IERR)
    elseif (tag(1:3) .eq. 'Q02') then
       !   Sample input 
       !   '+Q02    TiB Gradient   0lin 1lin/lin 2p/a   99file      '     0
       call divrd(CIOPTL,.TRUE., 0,.TRUE.,99,'TIB GRADIENT OPTION  ',IERR)
    elseif (tag(1:3) .eq. 'Q03') then
       !   Sample input 
       !   '+Q03 TN???? Te,Ti Flatten Option                         '     0
       call divrd(cflatopt,.true.,0,.true.,3,'Flatten Te,i option',ierr)
    elseif (tag(1:3) .eq. 'Q04') then
       !   Sample input 
       !   '+Q04 TN1447 Te flat upstream for S > SMAX * Teg cutoff = '    0.0
       call divrd(ctegcut,.FALSE.,0.0,.TRUE.,0.5,'TEB GRAD. CUTOFF',ierr)
    elseif (tag(1:3) .eq. 'Q05') then
       !   Sample input 
       !   '+Q05 TN1447 Ti flat upstream for S > SMAX * Tig cutoff = '    0.0
       call divrd(ctigcut,.FALSE.,0.0,.TRUE.,0.5,'TIB GRAD. CUTOFF',ierr)
    elseif (tag(1:3) .eq. 'Q06') then
       !   Sample input 
       !   '+Q06    Temperature of electrons at 0      TeB0  (eV)   '   20.0
       call divrd(CTEB0 ,.TRUE. ,0.0,.FALSE.,0.0,'TEMP AT 0   TEB0  ',IERR)
    elseif (tag(1:3) .eq. 'Q07') then
       !   Sample input 
       !   '+Q07    Temperature of electrons at plates TeBP  (eV)   '   20.0
       call divrd(CTEBP ,.TRUE. ,0.0,.FALSE.,0.0,'PLATES TEMP TEBP  ',IERR)
    elseif (tag(1:3) .eq. 'Q08') then
       !   Sample input 
       !   '+Q08    Temperature outer TeB step         TeBout(eV)   '    0.0
       call divrd(CTEBOU,.TRUE. ,0.0,.FALSE.,0.0,'TEMP STEP   TEBOUT',IERR)
    elseif (tag(1:3) .eq. 'Q09') then
       !   Sample input 
       !   '+Q09    Temperature inner TeB step         TeBin (eV)   '  100.0
       call divrd(CTEBIN,.TRUE. ,0.0,.FALSE.,0.0,'TEMP STEP   TEBIN ',IERR)
    elseif (tag(1:3) .eq. 'Q10') then
       !   Sample input 
       !   '+Q10    Temperature of trapped plasma      TeBt  (eV)   '   10.0
       call divrd(CTEBT ,.TRUE. ,0.0,.FALSE.,0.0,'TEMP TRAP   TEBT  ',IERR)
    elseif (tag(1:3) .eq. 'Q11') then
       !   Sample input 
       !   '+Q11 TN1278 Te exp decay step in trap      TeBouP       '    0.01
       call divrd(CTEBOUP,.TRUE.,0.0,.FALSE.,0.0,'TEMP DECAY TEBOUTP',IERR)
    elseif (tag(1:3) .eq. 'Q12') then
       !   Sample input 
       !   '+Q12    Temperature gradient factor        feBL1        '    0.0
       call divrd(CFEBL1,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FEBL1 ',IERR)
    elseif (tag(1:3) .eq. 'Q13') then
       !   Sample input 
       !   '+Q13    Temperature gradient factor        feBL2        '    0.0
       call divrd(CFEBL2,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FEBL2 ',IERR)
    elseif (tag(1:3) .eq. 'Q14') then
       !   Sample input 
       !   '+Q14    Temperature gradient factor        feBt         '    1.0
       call divrd(CFEBT ,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FEBT  ',IERR)
    elseif (tag(1:3) .eq. 'Q15') then
       !   Sample input 
       !   '+Q15    Temperature gradient factor        feB2         '    1.0
       call divrd(CFEB2 ,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FEB2  ',IERR)
    elseif (tag(1:3) .eq. 'Q16') then
       !   Sample input 
       !   '+Q16    Temperature of ions at 0           TiB0  (eV)   '   20.0
       call divrd(CTIB0 ,.TRUE. ,0.0,.FALSE.,0.0,'TEMP AT 0   TIB0  ',IERR)
    elseif (tag(1:3) .eq. 'Q17') then
       !   Sample input 
       !   '+Q17    Temperature of ions at plates      TiBP  (eV)   '   20.0
       call divrd(CTIBP ,.TRUE. ,0.0,.FALSE.,0.0,'PLATES TEMP TIBP  ',IERR)
    elseif (tag(1:3) .eq. 'Q18') then
       !   Sample input 
       !   '+Q18    Temperature outer TiB step         TiBout(eV)   '    0.05
       call divrd(CTIBOU,.TRUE. ,0.0,.FALSE.,0.0,'TEMP STEP   TIBOUT',IERR)
    elseif (tag(1:3) .eq. 'Q19') then
       !   Sample input 
       !   '+Q19    Temperature inner TiB step         TiBin (eV)   '  250.0
       call divrd(CTIBIN,.TRUE. ,0.0,.FALSE.,0.0,'TEMP STEP   TIBIN ',IERR)
    elseif (tag(1:3) .eq. 'Q20') then
       !   Sample input 
       !   '+Q20    Temperature of trapped plasma      TiBt  (eV)   '   10.0
       call divrd(CTIBT ,.TRUE. ,0.0,.FALSE.,0.0,'TEMP TRAP   TIBT  ',IERR)
    elseif (tag(1:3) .eq. 'Q21') then
       !   Sample input 
       !   '+Q21 TN1278 Ti exp decay step in trap      TiBouP       '    0.01
       call divrd(CTIBOUP,.TRUE.,0.0,.FALSE.,0.0,'TEMP DECAY TIBOUTP',IERR)
    elseif (tag(1:3) .eq. 'Q22') then
       !   Sample input 
       !   '+Q22    Temperature gradient factor        fiBL1        '    0.0
       call divrd(CFIBL1,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FIBL1 ',IERR)
    elseif (tag(1:3) .eq. 'Q23') then
       !   Sample input 
       !   '+Q23    Temperature gradient factor        fiBL2        '    0.0
       call divrd(CFIBL2,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FIBL2 ',IERR)
    elseif (tag(1:3) .eq. 'Q24') then
       !   Sample input 
       !   '+Q24    Temperature gradient factor        fiBt         '    1.0
       call divrd(CFIBT ,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FIBT  ',IERR)
    elseif (tag(1:3) .eq. 'Q25') then
       !   Sample input 
       !   '+Q25    Temperature gradient factor        fiB2         '    1.0
       call divrd(CFIB2 ,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FIB2  ',IERR)
    elseif (tag(1:3) .eq. 'Q26') then
       !   Sample input 
       !   '+Q26    Density at 0                       NB0   (m**-3)'    1.0E19
       call divrd(CNB0  ,.TRUE. ,0.0,.FALSE.,0.0,'DENS AT 0   NB0   ',IERR)
    elseif (tag(1:3) .eq. 'Q27') then
       !   Sample input 
       !   '+Q27 TN1264 Density at plates              NEBP  (m**-3)'    1.0e19
       call divrd(CNEBP ,.TRUE. ,0.0,.FALSE.,0.0,'PLATES DENS NEBP  ',IERR)
    elseif (tag(1:3) .eq. 'Q28') then
       !   Sample input 
       !   '+Q28    Density outer NB step              NBout (m**-3)'    0.03
       call divrd(CNBOUT,.TRUE. ,0.0,.FALSE.,0.0,'DENS STEP   NBOUT ',IERR)
    elseif (tag(1:3) .eq. 'Q29') then
       !   Sample input 
       !   '+Q29    Density inner NB step              NBin  (m**-3)'    5.0E18
       call divrd(CNBIN ,.TRUE. ,0.0,.FALSE.,0.0,'DENS STEP   NBIN  ',IERR)
    elseif (tag(1:3) .eq. 'Q30') then
       !   Sample input 
       !   '+Q30    Density of trapped plasma          NBt   (m**-3)'    1.0E19
       call divrd(CNBT  ,.TRUE. ,0.0,.FALSE.,0.0,'DENS TRAP   NBT   ',IERR)
    elseif (tag(1:3) .eq. 'Q31') then
       !   Sample input 
       !   '+Q31 TN1278 Ne exp decay step in trap      NBouP        '    0.01
       call divrd(CNBOUP,.TRUE. ,0.0,.FALSE.,0.0,'DENS DECAY  CNBOUP',IERR)
    elseif (tag(1:3) .eq. 'Q32') then
       !   Sample input 
       !   '+Q32 TN1347 Langmuir Probe Switch     0=Nb  1=Isat       '     1
       !     ENTER LANGMUIR PROBE DATA ALONG PLATES FOR TEMPERATURE GRADIENT
       !     OPTIONS 3 AND 4, AND PLASMA DECAY OPTIONS 2 AND 3. IN THE
       !     CASE OF THE PLASMA DECAY OPTIONS THE DATA ENTERED HERE WILL
       !     BE FILLED IN FOR ALL POINTS OF THE PLASMA. THUS GIVING A
       !     UNIFORM CONDITION UNLESS A TEMPERATURE GRADIENT OPTION
       !     MODIFIES THE TEMPERATURE.
       !
       !     FOR THOSE OPTIONS WITH ONLY ONE SET OF DATA THE INNER
       !     SPECIFICATION IS USED.
       !
       !     NOTE LPDATI :- LANGMUIR PROBE DATA INNER ...
       !
       !     THE LANGMUIR PROBE DATA SWITCH IDENTIFIES THE THIRD QUANTITY
       !     AS EITHER THE DENSITY OR THE Isat (PROBE SATURATION CURRENT)
       !
       !     THE DATA IS ORDERED AS  RING #,TEBP,TIBP,NBP (or ISAT)
       !
       call divrd(lpdatsw,.TRUE. , 0 ,.TRUE., 2 ,'LP DATA SWITCH', IERR)
    elseif (tag(1:3) .eq. 'Q33') then
       !   Sample input 
       !   '+Q33    OUTER Target Data Multipliers (Te,Ti,Nb):'  1.0   1.0   1.0
       call divrd(te_mult_i,ti_mult_i,n_mult_i,.TRUE.,0.0,.FALSE.,0.0,'Inner Multipliers       ',IERR)
    elseif (tag(1:3) .eq. 'Q34') then
       !   Sample input 
       !   '+Q34 ' 'Probe data at outer plate (opt4) or total (opt3)'
       !   '    Ring , TeBP , TiBP , NBP : Number of rows - '          0
       !
       !     Outer
       !
       call divrd(LPDATI,NLPDATI,MAXINS,-MACHHI,MACHHI,.FALSE.,0.0,MACHHI,3,'SET OF L. PROBE DATA INNER',IERR)
    elseif (tag(1:3) .eq. 'Q35') then
       !   Sample input 
       !   '+Q35    INNER Target Data Multipliers (Te,Ti,Nb):'  1.0   1.0   1.0
       call divrd(te_mult_o,ti_mult_o,n_mult_o,.TRUE.,0.0,.FALSE.,0.0,'Outer Multipliers       ',IERR)
    elseif (tag(1:3) .eq. 'Q36') then
       !   Sample input 
       !   '+Q36 ' 'Probe data at inner plate          (T grad opt4)'
       !   '    Ring , TeBP , TiBP , NBP : Number of rows - '          0
       call divrd(LPDATO,NLPDATO,MAXINS,-MACHHI,MACHHI,.FALSE.,0.0,MACHHI,3,'SET OF L. PROBE DATA OUTER',IERR)
    elseif (tag(1:3) .eq. 'Q37') then
       !   Sample input 
       !   '+Q37 ' 'CORE Plasma Data - for Core Options 1,2 and 3'
       !   '    Ring , TeB , TiB , NB , Vb : Number of rows - '       0
       !
       !     Read in data for core rings - data varies depending on
       !     core option selected.
       !
       call divrd(coredat,Ncoredat,MAXINS,-MACHHI,MACHHI,.FALSE.,-machhi,MACHHI,4,'SET OF CORE DATA',IERR)
    elseif (tag(1:3) .eq. 'Q38') then
       !   Sample input 
       !   '+Q38    Inboard plasma flow velocity       Vhin  (m/s)  '    0.0
       call divrd(CVHIN, .FALSE.,0.0,.FALSE.,0.0,'INB. PLASMA FLOW V',IERR)
    elseif (tag(1:3) .eq. 'Q39') then
       !   Sample input 
       !   '+Q39    Inboard electric field             Ein   (V/m)  '    0.0
       call divrd(CEIN , .FALSE.,0.0,.FALSE.,0.0,'INBOARD ELEC FIELD',IERR)
    elseif (tag(1:3) .eq. 'Q40') then
       !   Sample input 
       !   '+Q40    Outboard plasma flow vel  (SOL5)   Vhout (m/s)  '    0.0
       call divrd(CVHOUT,.FALSE.,0.0,.FALSE.,0.0,'OUT. PLASMA FLOW V',IERR)
    elseif (tag(1:3) .eq. 'Q41') then
       !   Sample input 
       !   '+Q41    Outboard electric field   (SOL5)   Eout  (V/m)  '    0.0
       call divrd(CEOUT ,.FALSE.,0.0,.FALSE.,0.0,'OUTBRD ELEC FIELD ',IERR)





    else if (tag(1:3) == 'Q42') then 
       ! ----------------------------------------------------------------------- 

       !     TAG Q42 : SHEATH TEMPERATURE - INNER JET/OUTER SONNET 


       !     Note: the tag line precedes a standard DIVIMP array input of 
       !           three lines. 

       !     - specifies a temperature value to be used instead of the 
       !       target value for the sheath energy calculations. 

       !     INPUT IS:  IR  TE 

       CALL RDRARN(sheath_vali,nsheath_vali,MAXNRS, &
            -MACHHI,MACHHI,.FALSE., &
            -machhi,MACHHI,1,'INNER/OUTER SHEATH POTENTIAL', &
            IERR) 



    else if (tag(1:3) == 'Q43') then 
       ! ----------------------------------------------------------------------- 

       !     TAG Q43 : SHEATH TEMPERATURE - OUTER JET/INNER SONNET 

       !     Note: the tag line precedes a standard DIVIMP array input of 
       !           three lines. 

       !     - specifies a temperature value to be used instead of the 
       !       target value for the sheath energy calculations. 

       !     INPUT IS:  IR  TE 

       CALL RDRARN(sheath_valo,nsheath_valo,MAXNRS, &
            -MACHHI,MACHHI,.FALSE., &
            -machhi,MACHHI,1,'OUTER/INNER SHEATH POTENTIAL', &
            IERR) 


    else if (tag(1:3) == 'Q44') then 
       ! ----------------------------------------------------------------------- 

       !     TAG Q44 : CORE PLASMA PROFILES AS A FUNCTION OF PSIN 


       !     Note: the tag line precedes a standard DIVIMP array input of 
       !           three lines. 

       !     - specifies the core plasma data as a function of PSIN 
       !       input line should be:  PSIN   Te    Ti    Ne    Vb 
       !       which will be linearly interpolated and overwrite data in Q37 if any 

       !     INPUT IS:  PSIN TE TI NE VB 


       CALL RDRARN_ALLOC(coreprofile,ncoreprofile, &
            -MACHHI,MACHHI,.FALSE., &
            -machhi,MACHHI,4,'CORE PLASMA PROFILES BY PSIN', &
            IERR) 

    else if (tag(1:3) == 'Q45') then 
       ! ----------------------------------------------------------------------- 

       !     TAG Q45 : DELTA PSI VALUE TO SHIFT THE INPUT CORE PROFILE 


       !     INPUT IS:  DELTA_PSIN_CORE 

       CALL ReadR(line,delta_psin_core,-machhi,machhi, &
            'PSIN SHIFT FOR CORE PROFILES') 



       !===================================================
       !
       ! TAG Series R
       !
       !===================================================
       !
       !  Inputs for TAG series R
       !
    elseif (tag(1:3) .eq. 'R01') then
       !   Sample input 
       !   '+R01      DETACHED PLASMA: Length Scaling Switch 0-S 1-P'    0
       call divrd(s21refsw,.true.,0,.true.,3,'S21 Length Ref Switch',ierr)
    elseif (tag(1:3) .eq. 'R02') then
       !   Sample input 
       !   '+R02 TN988 Detached Plasma Model: Te/Te0 at L1 O/I '  1.0    1.0
       call divrd(terat,terati,.true.,0.0,.false.,0.0,'Te Ratio at L1 O/I',ierr)
    elseif (tag(1:3) .eq. 'R03') then
       !   Sample input 
       !   '+R03 TN1439   Ti/Ti0 at Position L1            O/I '  1.0    1.0
       call divrd(tirat,tirati,.true.,0.0,.false.,0.0,'Ti Ratio at L1 O/I',ierr)
    elseif (tag(1:3) .eq. 'R04') then
       !   Sample input 
       !   '+R04 TN988    Ne/Ne0 at Position L1            O/I ' 10.0   10.0
       call divrd(nrat,nrati,.false.,0.0,.false.,0.0,'Ne Ratio at L1 O/I',ierr)
    elseif (tag(1:3) .eq. 'R05') then
       !   Sample input 
       !   '+R05 TN1496   Exponent for Ne Equation 1.0=lin O/I '  1.0    1.0
       call divrd(nalph,nalphi,.true.,0.0,.false.,0.0,'Ne Exponent at L1 O/I',ierr)
    elseif (tag(1:3) .eq. 'R06') then
       !   Sample input 
       !   '+R06 TN988    Qrad/Q0                          O/I '  5.0    5.0
       call divrd(qrat,qrati, .false.,0.0,.false.,0.0,'Emitted energy rat',ierr)
    elseif (tag(1:3) .eq. 'R07') then
       !   Sample input 
       !   '+R07 TN988    L1/SMAX ratio                    O/I '  0.1    0.1
       call divrd(l1rat,l1rati,.true.,0.0,.false.,0.0,'L1*SMAX boundary 1',ierr)
    elseif (tag(1:3) .eq. 'R08') then
       !   Sample input 
       !   '+R08 TN988    L2/SMAX ratio                    O/I '  0.2    0.2
       call divrd(l2rat,l2rati,.true.,0.0,.false.,0.0,'L2*SMAX boundary 2',ierr)
    elseif (tag(1:3) .eq. 'R09') then
       !   Sample input 
       !   '+R09 TN988    LV/SMAX ratio                    O/I '  0.2    0.2
       call divrd(lvrat,lvrati,.true.,0.0,.false.,0.0,'V->0 at lvrat*SMAX',ierr)
    elseif (tag(1:3) .eq. 'R10') then
       !   Sample input 
       !   '+R10 TN       Velocity multiplier at L1        O/I '  1.0    1.0
       call divrd(vbmult,vbmulti,.true.,0.0,.false.,0.0,'Vmult for Region A',ierr)
    elseif (tag(1:3) .eq. 'R11') then
       !   Sample input 
       !   '+R11  ' 'TN    DETACHED Plasma Specifications on a by ring basis'
       !   'TN    INNER - IR TER TIR NR QR L1R L2R LVR VBM - N= '    0
       call divrd(S21PARMI,NS21I,MAXNRS,-machhi,machhi,.FALSE.,-MACHHI,MACHHI,9,'SOL21 PARAMS INNER',IERR)
    elseif (tag(1:3) .eq. 'R12') then
       !   Sample input 
       !   '+R12 ' 'TN    DETACHED Plasma Specifications on a by ring basis'
       !   'TN    OUTER - IR TER TIR NR QR L1R L2R LVR VBM - N= '    0
       call divrd(S21PARMO,NS21O,MAXNRS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,9,'SOL21 PARAMS OUTER',IERR)

    else if (tag(1:3) == 'R13') then 
       ! ----------------------------------------------------------------------- 

       !     TAG R13 : SOL21 Additonal parameters - INNER JET/OUTER SONNET 

       !     Note: the tag line precedes a standard DIVIMP array input of 
       !           three lines. 

       !     - defines two additional points for linear fitting of N within 
       !       region A. 

       !     INPUT IS: IR L1A L1B NR1A NR1B TER1A TER1B TIR1A TIR1B 

       CALL RDRARN(aux_s21parmi,aux_ns21i,MAXNRS, &
            -MACHHI,MACHHI,.FALSE., &
            -machhi,MACHHI,8,'EXTRA SOL 21 INNER PARAMETERS', &
            IERR) 


    else if (tag(1:3) == 'R14') then 
       ! ----------------------------------------------------------------------- 

       !     TAG R14 : SOL21 Additonal parameters - OUTER JET/INNER SONNET 


       !     Note: the tag line precedes a standard DIVIMP array input of 
       !           three lines. 

       !     - defines two additional points for linear fitting of N within 
       !       region A. 

       !     INPUT IS: IR L1A L1B NR1A NR1B TER1A TER1B TIR1A TIR1B 

       CALL RDRARN(aux_s21parmo,aux_ns21o,MAXNRS, &
            -MACHHI,MACHHI,.FALSE., &
            -machhi,MACHHI,8,'EXTRA SOL 21 OUTER PARAMETERS', &
            IERR) 


    else if (tag(1:3) == 'R15') then 
       ! ----------------------------------------------------------------------- 

       !     TAG R15 : SOL21 extra velocity factor ... used to reduce flow 
       !               velocity near target. 

       CALL ReadR(line,aux_vel21,0.0,HI,'Modify near target velocity in SOL21') 




       !===================================================
       !
       ! TAG Series S
       !
       !===================================================
       !
       !  Inputs for TAG series S
       !
    elseif (tag(1:3) .eq. 'S01') then
       !   Sample input 
       !   '+S01    On-AXIS Toroidal B-field value                  '    2.0
       call divrd(cbphi,.true.,0.0,.false.,0.0,'ON-AXIS B-FIELD',    ierr)
    elseif (tag(1:3) .eq. 'S02') then
       !   Sample input 
       !   '+S02    Mass of plasma ions                Mb           '    2.0
       call divrd(CRMB, .TRUE. ,0.1 ,.FALSE.,0.0,'PLASMA ION MASS  ',IERR)
    elseif (tag(1:3) .eq. 'S03') then
       !   Sample input 
       !   '+S03    Charge on plasma ions              Zb           '    1
       call divrd(CIZB, .TRUE. , 1  ,.FALSE., 0 ,'PLASMA ION CHARGE',IERR) 
       RIZB = REAL (CIZB)
   elseif (tag(1:3) .eq. 'S04') then
       !   Sample input 
       !   '+S04    Mass of impurity ions              Mi           '   12.0
       call divrd(CRMI,  .TRUE. ,0.1,.FALSE.,0.0,'IMPURITY ION MASS', IERR)
    elseif (tag(1:3) .eq. 'S05') then
       !   Sample input 
       !   '+S05    Atomic number of impurity ions     Zi           '      6
       call divrd(CION,  .TRUE. , 1 ,.FALSE., 0 ,'IMP ATOMIC NUMBER', IERR)
    elseif (tag(1:3) .eq. 'S06') then
       !   Sample input 
       !   '+S06    Initial temperature                Tem1  (eV)   '    0.5
       call divrd(CTEM1, .FALSE.,0.0,.FALSE.,0.0,'INITIAL TEMP   ',   IERR)
    elseif (tag(1:3) .eq. 'S07') then
       !   Sample input 
       !   '+S07    Initial temperature (2)            Tem2  (eV)   '    0.0
       call divrd(CTEM2 ,.FALSE.,0.0,.FALSE.,0.0,'INITIAL TEMP (2)',  IERR)
    elseif (tag(1:3) .eq. 'S08') then
       !   Sample input 
       !   '+S08    Initial R position of impurity     R0    (m)    '    1.4
       call divrd(CXSC,  .FALSE.,0.0,.FALSE.,0.0,'INITIAL R POSITION',IERR)
    elseif (tag(1:3) .eq. 'S09') then
       !   Sample input 
       !   '+S09    Initial Z position of impurity     Z0    (m)    '    0.9
       call divrd(CYSC,  .FALSE.,0.0,.FALSE.,0.0,'INITIAL Z POSITION',IERR)
    elseif (tag(1:3) .eq. 'S10') then
       !   Sample input 
       !   '+S10    Operation Mode  1 Time-Dependent  2 Steady-State'      2
       !     jdemod - remove upper bounds checks for maxizs and maximp due
       !              to dynamical allocation
       call divrd(IMODE, .TRUE.,  0, .TRUE.,   2  ,'OPERATION MODE',  IERR)
    elseif (tag(1:3) .eq. 'S11') then
       !   Sample input 
       !   '+S11    Number of impurity ions to be followed          '    100
       call divrd(NIMPS, .TRUE.,  1, .FALSE.,MAXIMP,'NO OF IONS',     IERR)
    elseif (tag(1:3) .eq. 'S12') then
       !   Sample input 
       !   '+S12 TN487 Number of Supplementary Neutrals to Launch   '      0
       call divrd(NIMPS2,.TRUE.,0,.FALSE.,MAXIMP-NIMPS,'NUM SUP IONS',IERR)
    elseif (tag(1:3) .eq. 'S13') then
       !   Sample input 
       !   '+S13    Quantum iteration time for atoms   fsrate (s)   '    1.0E-8
       call divrd(FSRATE,.TRUE.,1.E-10,.TRUE.,1.0,'NEUT ITERATE TIME',IERR)
    elseif (tag(1:3) .eq. 'S14') then
       !   Sample input 
       !   '+S14    Quantum iteration time for ions    qtim   (s)   '    1.0E-8
       call divrd(QTIM,  .TRUE.,1.E-10,.TRUE.,1.0,'DIV ITERATE TIME', IERR)
    elseif (tag(1:3) .eq. 'S15') then
       !   Sample input 
       !   '+S15 T   CPU time limit in seconds          cpulim (s)   '  80000.0
       !---- READ IN TIME DEPENDENT DATA  (NOTE 128)
       !
       call divrd(CPULIM,.TRUE.,0.0,.FALSE., 0.0 ,'CPU TIME LIMIT'  , IERR)
    elseif (tag(1:3) .eq. 'S16') then
       !   Sample input 
       !   '+S16 ' 'Average Dwell Times (s) for each state 0,1,2..'
       !   '        Number of dwell times given below :-'    0
       call divrd(DWELTS,NQS,MAXIZS+1,0.0,machhi,.FALSE.,'DWELT',IERR)
       !call divrd(DWELTS(0),NQS,MAXIZS+1,0.0,machhi,.FALSE.,'DWELT',IERR)
       ! jdemod:
       ! DWELTS has been allocated to a size of -1:maxizs+1
       ! The code originally loaded this array from element 0 but the revised code starts with the first array element which is -1
       ! This means that the stored values will have been assigned starting at element -1 but are expected to start at index 0.
       ! This is corrected here by moving up the elements in the array by one
       ! DWELTS(-1) is supposed to equal DWELTS(0) so this element is not changed.
       ! Move DWELTS entries up by one array element
       do in = nqs,0
          dwelts(in) = dwelts(in-1)
       end do
    elseif (tag(1:3) .eq. 'S17') then
       !   Sample input 
       !   '+S17 ' 'Dwell Time Factors for time-dependent analysis'
       !   '        Number of dwell time factors :-'  0
       call divrd(DWELFS,NTS,MAXNTS,  0.0,machhi,.TRUE.,'T FACTORS',IERR)
    elseif (tag(1:3) .eq. 'S18') then
       !   Sample input 
       !   '+S18 T   Maximum dwell time for steady state (s)         '    0.5
       call divrd(ctimmax,.true.,0.0,.false.,10.0,'DWELL TIME LIMIT ',ierr)
    elseif (tag(1:3) .eq. 'S19') then
       !   Sample input 
       !   '+S19    Random number seed  (0 generate new seed)       '      0
       call divrd(CISEED,.FALSE., 0 ,.FALSE., 0 ,'RANDOM NUMBER SEED',IERR)
    elseif (tag(1:3) .eq. 'S20') then
       !   Sample input 
       !   '+S20    Number of Iterations                            '      1
       call divrd(NITERS,.TRUE. , 1 ,.TRUE. ,20 ,'NO. OF ITERATIONS ',IERR)
    elseif (tag(1:3) .eq. 'S21') then
       !   Sample input 
       !   '+S21    SOLTEST - 0.0 run normally -1.0 test SOL opt    '     0.0
       call divrd(ctestsol,.true.,-1.0,.false.,1.0,'TEST SOL OPT ONLY',ierr)

    else if (tag(1:3) == 'S22') then 
       ! ----------------------------------------------------------------------- 

       !     TAG S22 : Option to turn on V/A flag debugging in NEUT 
       !               OFF =0 = default 
       !               ON  >0 

       CALL ReadI(line,debug_neutv,0,1,'Neutral velocity debugging') 

    else if (tag(1:3) == 'S23') then 
       ! ----------------------------------------------------------------------- 

       !     TAG S23 : Maximum energy to be used for binning the 
       !               velocity/energy in the debugging option code 

       CALL ReadR(line,debug_neutv_einmax,0.0,HI, &
            'Max EIN for NEUT V/A binning for debugging') 

    else if (tag(1:3) == 'S24') then 
       ! ----------------------------------------------------------------------- 

       !     TAG S24 : Number of bins to be used in Velocity debugging code 
       !               Maximum of 20000. 

       CALL ReadI(line,debug_neutv_nbins,1,20000, &
            'Neutral velocity debugging - number of bins') 

       !...    Using 'S' prefix: 
    else if (tag(1:3) == 'S70') then 
       CALL ReadR(line,osmbulkn,-HI,HI,'bluk n over-ride') 
    else if (tag(1:3) == 'S71') then 
       CALL ReadR(line,osmbulkte,-HI,HI,'bluk electron T  over-ride') 
    else if (tag(1:3) == 'S72') then 
       CALL ReadR(line,osmbulkti,-HI,HI,'bluk ion T over-ride') 
    else if (tag(1:3) == 'S73') then 
       CALL ReadR(line,osmbulkv,-HI,HI,'bluk plasma v over-ride') 
    else if (tag(1:3) == 'S74') then 
       !... 
       READ(line(7:9),*) s28mode 
       if (s28mode == 1.0.or.s28mode == 1.1) then 
          CALL RDRARN(osms28,osmns28,MAXNRS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,7,'SOL28 specifications 1.x',IERR) 
          DO i1 = 1, osmns28 
             osms28(i1,9)  = 0.0 
             osms28(i1,10) = 0.0 
             WRITE(fp,'(A,1P,1P,8E10.2,0P)') 'S28:',(osms28(i1,i2),i2 = 1,8) 
          end do
       else if (s28mode == 2.0) then 
          CALL RDRARN(osms28,osmns28,MAXNRS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,7,'SOL28 specifications 2.0',IERR) 
          DO i1 = 1, osmns28 
             osms28(i1,9)  = 0.0 
             osms28(i1,10) = 0.0 
             WRITE(fp,'(A,1P,1P,8E10.2,0P)') 'S28:',(osms28(i1,i2),i2 = 1,8) 
          end do
       else if (s28mode == 2.1) then 
          !...      Added ability to limit application of SOL28 parameters to specific 
          !         rings: 
          CALL RDRARN(osms28,osmns28,MAXNRS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,9,'SOL28 specifications 2.1',IERR) 
          DO i1 = 1, osmns28 
             WRITE(fp,'(A,1P,1P,10E10.2,0P)') &
                  'S28:',(osms28(i1,i2),i2 = 1,10) 
          end do
       else if (s28mode >= 3.0) then 
          !...      Improved: 
          CALL RDRARN(osms28,osmns28,MAXNRS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,11,'SOL28 specifications 3.0',IERR) 
          DO i1 = 1, osmns28 
             WRITE(fp,'(A,1P,1P,12E10.2,0P)') &
                  'S28:',(osms28(i1,i2),i2 = 1,12) 
          end do
       ELSE 
          CALL ER('RUI','Unsupported version for S74 OSMS28',*99) 
       end if
    else if (tag(1:3) == 'S75') then 
       !... 
       CALL Read2I(line,s28ion,s28ionpfz,-5,20,'SOL28 ion src option') 
    else if (tag(1:3) == 'S83') then 
       !... 
       CALL ReadR(line,s28ionbound,-99.0,5.0,'SOL28 ion limit') 
    else if (tag(1:3) == 'S84') then 
       !... 
       CALL Read2I(line,s28cfp,s28cfppfz,0,25,'SOL28 cross-f ion opt') 
    else if (tag(1:3) == 'S76') then 
       !... 
       CALL Read2I(line,s28mom,s28mompfz,-5,20,'SOL28 momentum option') 
    else if (tag(1:3) == 'S85') then 
       !... 
       CALL Read2I(line,s28cfm,s28cfmpfz,0,1,'SOL28 cross-f mom opt') 
    else if (tag(1:3) == 'S98') then 
       !... 
       CALL Read2I(line,s28te,s28tepfz,0,21,'SOL28 Te option') 
    else if (tag(1:3) == 'S86') then 
       !... 
       CALL Read2I(line,s28ti,s28tipfz,0,20,'SOL28 Ti option') 
    else if (tag(1:3) == 'S87') then 
       !... 
       CALL ReadR(line,s28tiratio,0.1,10.0,'SOL28 Ti/Te ratio') 
    else if (tag(1:3) == 'S88') then 
       !... 
       CALL Read2I(line,s28superdet,s28superdetpfz,-5,5,'SOL28 M opt') 
    else if (tag(1:3) == 'S77') then 
       !... 
       CALL ReadR(line,s28momfr,0.0,HI,'SOL28 momentum balance') 
    else if (tag(1:3) == 'S78') then 
       !... 
       CALL ReadI(line,s22pfzion,0,1,'feeding PFZ ionisation') 
    else if (tag(1:3) == 'S79') then 
       !... 
       CALL ReadR(line,setmach,0.0,1.0,'setting mach number in SOL24') 
    else if (tag(1:3) == 'S80') then 
       !... 
       CALL ReadR(line,muldensity,-HI,HI,'SOL24,28 density multi') 
    else if (tag(1:3) == 'S81') then 
       !... 
       CALL ReadR(line,addtemp,-10.0,99.0,'SOL24,28 density temp add') 
    else if (tag(1:3) == 'S82') then 
       !... 
       CALL ReadR(line,s28b34,0.0,1.0,'SOL28 Te profile exponent') 
    else if (tag(1:3) == 'S89') then 
       !... 
       CALL ReadR(line,s28b34det,0.0,1.0,'SOL28 det Te profile exp') 
    else if (tag(1:3) == 'S90') then 
       !... 
       CALL Read2I(line,s28probe,s28probepfz,0,999,'SOL28 probe data') 
    else if (tag(1:3) == 'S91') then 
       !... 
       CALL Read2I(line,s28rec,s28recpfz,0,20,'SOL28 recomb. opt') 
    else if (tag(1:3) == 'S92') then 
       !... 
       CALL ReadR(line,s28momTe,0.0,30.0,'SOL28 mom Te cutoff') 
    else if (tag(1:3) == 'S93') then 
       !... 
       READ(line(7:9),*) version 
       if (version == 1.0) then 
          CALL RDRARN(osms28over,osmns28over,MAXNRS,-MACHHI,MACHHI, &
               .FALSE.,-MACHHI,MACHHI,8,'SOL28 over-ride 1.0', &
               IERR) 
          DO i1 = 1, osmns28over 
             WRITE(fp,'(A,1P,9E10.2)') 'OVER:',(osms28over(i1,i2),i2 = 1,9) 
          end do
       ELSE 
          CALL ER('RUI','Unsupported version for S93 OSMS28OVER',*99) 
       end if
    else if (tag(1:3) == 'S94') then 
       !... 
       CALL ReadI(line,s28sym,0,1,'Symmeterize non-SOL28 rings') 
    else if (tag(1:3) == 'S95') then 
       !... 
       CALL Read2I(line,s28cfpnt,s28cfppfznt,0,7,'Near target cf src') 
    else if (tag(1:3) == 'S96') then 
       !... 
       CALL Read2I(line,s28cfpdrft,s28cfppfzdrft,-3,6,'Rad. drift') 
    else if (tag(1:3) == 'S97') then 
       !... 
       CALL Read2I(line,s28nemode,s28nemodepfz,-3,2,'upstream ne/jsat') 



       !===================================================
       !
       ! TAG Series T
       !
       !===================================================
       !
       !  Inputs for TAG series T
       !
    elseif (tag(1:3) .eq. 'T01') then
       !   Sample input 
       !   '+T01    Ionization     0/1/2old 3ei 4no 5ei/dis 6no/dis '     3
       call divrd(CIOPTA,.TRUE., 0,.TRUE., 6,'IONIZATION OPT       ',IERR)
    elseif (tag(1:3) .eq. 'T02') then
       !   Sample input 
       !   '+T02    Collision opt  0std 1inf 2zeff 3tpara           '    13
       !     All other Reiser input is unstructured and optional
       call divrd(CIOPTB,.TRUE., 0,.TRUE.,14,'COLLISION OPT        ',IERR)
    elseif (tag(1:3) .eq. 'T03') then
       !   Sample input 
       !   '+T03    REISER                                          '     0
       call divrd(CIOPTR,.TRUE., 0,.TRUE., 2,'REISER COLLISION OPT ',IERR)
    elseif (tag(1:3) .eq. 'T04') then
       !   Sample input 
       !   '+T04    Friction opt   0std 1inf 2tpara                 '     0
       call divrd(CIOPTC,.TRUE., 0,.TRUE., 4,'FRICTION OPT         ',IERR)
    elseif (tag(1:3) .eq. 'T05') then
       !   Sample input 
       !   '+T05    Heating opt    0std 1inf 2zero 3Ti              '     0
       call divrd(CIOPTD,.TRUE., 0,.TRUE., 3,'HEATING OPT          ',IERR)
    elseif (tag(1:3) .eq. 'T06') then
       !   Sample input 
       !   '+T06    CX Recomb opt  0off 1on 2Vcx                    '     0
       call divrd(CIOPTI,.TRUE., 0,.TRUE., 9,'CX RECOMB OPT        ',IERR)
    elseif (tag(1:3) .eq. 'T07') then
       !   Sample input 
       !   '+T07    Dperp option   0const 1vary                     '     0
       call divrd(CIOPTJ,.TRUE., 0,.TRUE., 4,'DPERP OPTION         ',IERR)
    elseif (tag(1:3) .eq. 'T08') then
       !   Sample input 
       !   '+T08 TN1272 Perpendicular step option 0-normal 1-core    '     3
       call divrd(cdiffopt,.true.,0,.true.,3,'PERP STEP PROB OPTION',ierr)
    elseif (tag(1:3) .eq. 'T09') then
       !   Sample input 
       !   '+T09 TN14?? Pinch Velocity Option 0=off 1=all 2=main SOL '     0
       call divrd(pinchopt,.true.,0,.true.,15,'PINCH VELOCITY OPTION',ierr)
    elseif (tag(1:3) .eq. 'T10') then
       !   Sample input 
       !   '+T10    TeB Grad Coeff 0off 1on                         '     1
       call divrd(CIOPTM,.TRUE., 0,.TRUE., 3,'TEB GRAD COEFF OPTION',IERR)
    elseif (tag(1:3) .eq. 'T11') then
       !   Sample input 
       !   '+T11    TiB Grad Coeff 0off 1on                         '     1
       call divrd(CIOPTN,.TRUE., 0,.TRUE., 3,'TIB GRAD COEFF OPTION',IERR)
    elseif (tag(1:3) .eq. 'T12') then
       !   Sample input 
       !   '+T12    T-GRAD Force Modification Function 0-off 1-UEDGE'     0
       call divrd(fgradopt,.true.,0,.true.,4,'GRAD FORCE MOD OPTION',ierr)
    elseif (tag(1:3) .eq. 'T13') then
       !   Sample input 
       !   '+T13 TN505 Poloidal Velocity Drift Option 0-off 1-on     '     0
       call divrd(CPDRFT,.TRUE., 0,.TRUE., 3,'POL. DRIFT OPTION    ',IERR)
    elseif (tag(1:3) .eq. 'T14') then
       !   Sample input 
       !   '+T14    Cross Field Diffusion factor       Dperp (m*m/s)'    0.3
       call divrd(CDPERP,.TRUE.,0.0,.FALSE.,0.0,'X DIFF DPERP       ',IERR)
    elseif (tag(1:3) .eq. 'T15') then
       !   Sample input 
       !   '+T15    Trap Cross Field Diffusion factor  Dperpt(m*m/s)'   -1.0
       !   if (cdperpt.le.0.0) cdperpt = cdperp
       ! move after input file read
       call divrd(CDPERPT,.FALSE.,0.0,.FALSE.,0.0,'X DIFF DPERP-TRAP',IERR)
    elseif (tag(1:3) .eq. 'T16') then
       !   Sample input 
       !   '+T16    Perpendicular Pinch Velocity        Vpinch (m/s)'    0.0
       call divrd(cVPINCH,.FALSE.,0.0,.FALSE.,0.0,'PERP PINCH VEL.  ',IERR)
    elseif (tag(1:3) .eq. 'T17') then
       !   Sample input 
       !   '+T17 TN505 Poloidal Drift Velocity (m/s)                 '   0.0
       call divrd(CDRFTV,.FALSE.,0.0,.FALSE.,0.0,'POL. DRIFT VEL.  ',IERR)
    elseif (tag(1:3) .eq. 'T18') then
       !   Sample input 
       !   '+T18 TN    Poloidal Drift Range F1*SMAX < S < F2*SMAX'    0.2 0.8
       call divrd(cdrftv_start,cdrftv_end,.false.,0.0,.false.,1.0,'Start and End of DRFTV',ierr)

    else if (tag(1:3) == 'T19') then 
       ! ----------------------------------------------------------------------- 

       !     TAGS T19 to T27 

       !     The following are all variables related to the REISER 
       !     force implementation. T10 to T16 allow for detailed control of the 
       !     active forces while T17 is the coulomb parameter and T18 is 
       !     a debug option to linearize the gradients at the midplane in the 
       !     case of SOL option 7. 


       CALL ReadI(line,aswitch,0,1,'Master Switch for force switches') 
    else if (tag(1:3) == 'T20') then 
       CALL ReadI(line,sk11,0,1,'FF switch') 
    else if (tag(1:3) == 'T21') then 
       CALL ReadI(line,sk12,0,1,'FIG switch') 
    else if (tag(1:3) == 'T22') then 
       CALL ReadI(line,sk13,0,1,'FVG switch') 
    else if (tag(1:3) == 'T23') then 
       CALL ReadI(line,sd11,0,1,'Maxwellian Velocity diffusion switch') 
    else if (tag(1:3) == 'T24') then 
       CALL ReadI(line,sd12,0,1,'TiGradient Velocity diffusion switch') 
    else if (tag(1:3) == 'T25') then 
       CALL ReadI(line,sd13,0,1,'V-Gradient Velocity diffusion switch') 
    else if (tag(1:3) == 'T26') then 
       CALL ReadR(line,coulomb_log,0.0,HI,'Coulomb Logarithm Value') 
    else if (tag(1:3) == 'T27') then 
       CALL ReadI(line,linearpeak,0,1,'Linear Peak Debug switch') 



    else if (tag(1:3) == 'T28') then 
       ! ----------------------------------------------------------------------- 

       !     TAG T28 

       !     T28 - pinch_loc_opt - this option specifies the region of the 
       !                           grid where the radial velocity in 
       !                           pinchopt 4 will be applied. 

       !         = 0 = Entire grid excluding PFZ (default) 
       !         = 1 = main SOL only 
       !         = 2 = Entire grid including PFZ 
       !               - this requires a sign change to the value assigned to 
       !                 Vr (or Vpinch depending on terminology) 
       !         = 3 = main SOL past X-point region ONLY 
       !         = 4 = main SOL past X-point + core 

       CALL ReadI(line,pinch_loc_opt,0,4,'Grid region for pinch') 

    else if (tag(1:3) == 'T29') then 
       ! ----------------------------------------------------------------------- 

       !     TAG T29 

       !     T29 - pinch_npdf, pinch_pdf - this option loads the probability 
       !           distribution function to be used when randomly 
       !           determining the value of the pinch/radial velocity at 
       !           each time step. 

       !     Note: the tag line precedes a standard DIVIMP array input of 
       !           three lines. 


       CALL RDRARN(pinch_pdf,pinch_npdf,MAXPTS, &
            -MACHHI,MACHHI,.TRUE., &
            -machhi,MACHHI,1,'RADIAL/PINCH V - PDF', &
            ierr) 

       !       Calculate the integral value of the PDF and assign 
       !       this to the normalization value. 

       !       Note: the integration assumes that P=0 for V<Vmin-1/2*dVbin 
       !             and V>Vmax+1/2*dVbin 
       !             so the implicit integration range is 
       !             [Vmin-1/2dV,Vmax+1/2dV] of the 
       !             input pdf values. 

       !       The PDF is input in the form: 

       !           Vel    PDF-value 

       !       and the values are linearly interpolated between the points on 
       !       the function to generate intermediate PRPBABILITIES. 

       !       Calculate the integral and store it at each point in the data 
       !       - take the total and assign it to pdf_norm_val 

       !       Set maximum and minimum velocities 

       npdf_data = pinch_npdf+2 

       pinch_pdf_data(1,1) = pinch_pdf(1,1) - &
            (pinch_pdf(2,1)-pinch_pdf(1,1))/2.0 
       pinch_pdf_data(1,2) = 0.0 

       pinch_pdf_data(npdf_data,1) = pinch_pdf(pinch_npdf,1)+ &
            (pinch_pdf(pinch_npdf,1)-pinch_pdf(pinch_npdf-1,1))/2.0 
       pinch_pdf_data(npdf_data,2) = 0.0 

       !       Copy probability information to data array. 

       do in = 1,pinch_npdf 
          pinch_pdf_data(in+1,1) = pinch_pdf(in,1) 
          pinch_pdf_data(in+1,2) = pinch_pdf(in,2) 
       end do

       !       Calculate integral 

       !       Integration starts at zero 

       pdf_norm_val = 0.0 
       pinch_pdf_data(1,3) = pdf_norm_val 

       do in = 1,npdf_data-1 

          !          Define the velocity bin sizes based on input 

          !          Calculate the integral by assuming the probability 
          !          at each point follows the specified input rather 
          !          than being constant for each dVbin. 

          !           if (in.eq.1) then 
          !              deltav2 = (pinch_pdf(in+1,1)-pinch_pdf(in,1))/2.0 
          !              deltav1 = deltav2 
          !           elseif (in.eq.pinch_npdf) then 
          !              deltav1 = (pinch_pdf(in,1)-pinch_pdf(in-1,1))/2.0 
          !              deltav2 = deltav1 
          !           else 
          !              deltav1 = (pinch_pdf(in,1)-pinch_pdf(in-1,1))/2.0 
          !              deltav2 = (pinch_pdf(in+1,1)-pinch_pdf(in,1))/2.0 
          !           endif 

          deltav1 = pinch_pdf_data(in+1,1) - pinch_pdf_data(in,1) 

          pdf_norm_val = pdf_norm_val &
               + (pinch_pdf_data(in+1,2)+pinch_pdf_data(in,2))/2.0 &
               * deltav1 

          !          Accumulate integral at the bin centers 
          !          Integrated probability is zero in first cell. 

          pinch_pdf_data(in+1,3) = pdf_norm_val 

       end do

       !       Normalize the integral of the PDF - used in random number selection 

       if (pdf_norm_val /= 0.0) then 

          do in = 1,npdf_data 

             pinch_pdf_data(in,3) = pinch_pdf_data(in,3)/pdf_norm_val 

             write(6,'(a,i4,4(1x,f14.8))') 'PINCH INT:',in, &
                  pinch_pdf_data(in,1),pinch_pdf_data(in,2), &
                  pinch_pdf_data(in,3), pdf_norm_val 

          end do

       end if

       !      Print out vr_pdf_int 

       !       do in = -40,50 
       !          vtest = 10.0 * in 
       !          res   = vr_pdf_int(vtest,-1) 
       !          write(6,'(a,i4,3(1x,f14.8))') 'PDF_INT:',in, 
       !     >       vtest,res,res/pdf_norm_val 
       !       end do 


    else if (tag(1:3) == 'T30') then 
       ! ----------------------------------------------------------------------- 

       !     TAG T30 

       !     T30 - pinch correlation time. A new radial velocity value will be 
       !           chosen periodically based on the value of this quantity. 
       !           The default value of 0.0 will result in a new velocity 
       !           being selected every time step. This value is specified 
       !           in second. On JET it is typically 5 to 20 microseconds. 


       CALL ReadR(line,pinch_correlation_time,0.0,HI, &
            'Pinch Correlation Time') 



    else if (tag(1:3) == 'T31') then 
       ! ----------------------------------------------------------------------- 

       !     TAG T31 

       !     T31 - Drift region - specifies the region to which poloidal 
       !           drifts should be applied. 
       !           1 - SOL + PFZ 
       !           2 - SOL only 
       !           3 - PFZ only 
       !           4 - CORE only 

       !           Other options can easily be added as needed - the default is 
       !           option 1. 


       CALL ReadI(line,drft_region,0,4,'Region over which to apply'// &
            'poloidal drifts') 


    else if (tag(1:3) == 'T32') then 
       ! ----------------------------------------------------------------------- 

       !     TAG T32 

       !     T32 - Drift Mach Option - Detailed drift velocity input on a ring 
       !           ring basis is specified as a mach number to be multiplied 
       !           by the sound speed at the top of the torus for each ring 

       !           Option 0 : OFF 
       !           Option 1 : CS calculated from 2*Te 
       !           Option 2 : CS calculated from Te+Ti 

       !           Default value is 0 - OFF - data is specified in terms of 
       !                                      velocity 

       CALL ReadI(line,drftvel_machopt,0,2,'Drift velocity specified'// &
            'by mach values') 

    else if (tag(1:3) == 'T33') then 
       ! ----------------------------------------------------------------------- 

       !     TAG T33 

       !     T33 - Detailed specifications of data ring by ring - this is 
       !           an array listing 
       !                      ring number       velocity/mach 
       !           Data does not need to be specified for each ring - the 
       !           default value will be applied instead. 

       !           The number of array elements is initialized to zero 


       !        Read in array data 

       CALL RDRARN(ringdrftvel,ndrftvel,MAXNRS, &
            real(1),real(maxnrs),.TRUE., &
            -machhi,MACHHI,1,'DETAILED RING DRIFT VELOCITY', &
            ierr) 



    else if (tag(1:3) == 'T34') then 
       ! ----------------------------------------------------------------------- 

       !     TAG T34 

       !     T34 - S displacement in 2D resulting from a perpendicular step in 
       !           paramagnetic "Z" direction in 3D - this actually moves the 
       !           particle onto an adjacent flux tube - however, since 
       !           DIVIMP is 2D - the effect is to actually move the particle 
       !           onto an adjacent identical flux tube at a different S location 
       !           - thus effectively giving a net S displacement. 

       !           The first approximation to this is to use 

       !           ds = cross_step * Btor/Bpol 

       !           A value of 0 for this option is OFF 
       !                      1 is ON 
       CALL ReadI(line,dperpz_opt,0,1,'3D Dperp Delts S Option') 



    else if (tag(1:3) == 'T35') then 
       ! ----------------------------------------------------------------------- 

       !     TAG T35 - related to poloidal drift options - T31,T32,T33 

       !     T35 - Drift region calculation option 

       !     Option 0: Input values are specified in terms of S 
       !     Option 1: Input values are specified in terms of P (poloidal distance) 
       !     Option 2: Input values are specified in terms of Z (single null only) 

       !     Default value is S para (option 0) 

       !     jdemod - I don't really understand commenting this out 
       !            - it defaults to a value of 0 so this will only stop the code 
       !              if the T35 line is in the input file 
       !            - other than that I think it works as intended so I am uncommenting 
       !              it .. its only purpose is to change the interpretation of other 
       !              inputs 

       ! slmod begin - *** TEMP *** 
       CALL ReadI(line,drft_distopt,0,2,'Drift velocity range'// &
            ' specification (S,P or Z)') 

       !        STOP 'OPTION TURNED OFF FOR NOW...' 
       ! slmod end 

       !        write(0,*) 'READIN: drft_distopt:',drft_distopt 



    else if (tag(1:3) == 'T36') then 
       ! ----------------------------------------------------------------------- 

       !     TAG T36 to T39 - options related to the implementation of 
       !                      impurity exb drifts 

       !     TAG T36 - potopt 

       !             - This option is used to determine the method of calculatng 
       !               plasma potential 
       !             - potopt = 0    Use 3xTe(0) at each target as the floating 
       !                              potential start start point 
       !             - potopt = 1    Import LP data listing the measured floating 
       !                              potential ... if imported data not available it 
       !                              defaults to option 0.(not yet implemented) 

       CALL ReadI(line,potopt,0,0,'Option for calculating the floating' &
            //' potential') 


    else if (tag(1:3) == 'T37') then 
       !     TAG 37 

       !     exb_rad_opt = 0 ... no exb radial drift is applied 
       !                 = 1 ... exb radial drift is turned on 

       CALL ReadI(line,exb_rad_opt,0,1,'ExB radial drift option') 


    else if (tag(1:3) == 'T38') then 
       !     TAG 38 

       !     exb_pol_opt = 0 ... no exb poloidal drift is applied 
       !                 = 1 ... exb poloidal drift is turned on 

       CALL ReadI(line,exb_pol_opt,0,1,'ExB poloidal drift option') 

    else if (tag(1:3) == 'T39') then 

       !     TAG 39 

       !     exb_scale - real number 
       !               - the basic function of this is to switch the sign of 
       !                 the ExB drift for cases of forward (+1.0) and reverse (-1.0) 
       !                 B-field orientation. However, it can also be used as a scaling 
       !                 factor if the drifts are found to be either too large or too 
       !                 small. 
       CALL ReadR(line,exb_scale,-HI,HI, &
            'ExB scaling factor (usually +/-1.0)') 



    else if (tag(1:3) == 'T40') then 

       ! ----------------------------------------------------------------------- 

       !     Force scaling factors - all default to 1.0 
       !     T40 to T44 
       !     T40 = Friction force scaling factor (SF_FRIC) 
       !     T41 = Ion temperature force scaling factor (SF_TI) 
       !     T42 = Electron temperature force scaling factor (SF_TE) 
       !     T43 = Electric field force scaling factor (SF_EF) 
       !     T44 = Velocity diffusion scaling factor (SF_VDIFF) 
       !     T45 = Scaling factor for TAU (sf_tau) 

       !     Defaults: 
       !     sf_fric = 1.0 
       !     sf_ti   = 1.0 
       !     sf_te   = 1.0 
       !     sf_ef   = 1.0 
       !     sf_vdiff= 1.0 
       !     sf_tau  = 1.0 

       !----------------------------------------------------------------------- 

       CALL ReadR(line,sf_fric,-HI,HI, &
            'Friction force scaling factor (def = 1.0)') 
    else if (tag(1:3) == 'T41') then 
       CALL ReadR(line,sf_ti,-HI,HI, &
            'Ti gradient force scaling factor (def = 1.0)') 
    else if (tag(1:3) == 'T42') then 
       CALL ReadR(line,sf_te,-HI,HI, &
            'Te gradient force scaling factor (def = 1.0)') 
    else if (tag(1:3) == 'T43') then 
       CALL ReadR(line,sf_ef,-HI,HI, &
            'Electric field force scaling factor (def = 1.0)') 
    else if (tag(1:3) == 'T44') then 
       CALL ReadR(line,sf_vdiff,-HI,HI, &
            'Velocity diffusion scaling factor (def = 1.0)') 
    else if (tag(1:3) == 'T45') then 
       CALL ReadR(line,sf_tau,-HI,HI, &
            'TAU scaling factor (def = 1.0)') 


    else if (tag(1:3) == 'T46') then 
       ! ----------------------------------------------------------------------- 

       !    T46 Velocity based temperature calculation option 

       !    0  = default 
       !    1+ = other options 

       CALL ReadI(line,ti_calc_opt,0,3,'Impurity Ti Calculation Opt') 

    else if (tag(1:3) == 'T47') then 
       ! ----------------------------------------------------------------------- 

       !    T47 Coulomb logarithm calculation options 

       !     0  = default = constant (default value = 15.0) 
       !     1  = 30.0 - 0.5 * LOG(ni) + 1.5 * LOG(ti)  [HC code - ] 
       !          Originally in Sivukhin, D.V., Coulomb collisions in a fully ionized plasma in 
       !          Review of Plasma Physics (Consultation Bureau, New York, 1966) Vol. 4, p.88. 

       !     2  = 17.3 - 0.5*LOG(n/1.0E20) + 1.5*LOG(t/1000.0)  [LIM code] 
       !     3  = log(1.5e13 * t**(1.5) / sqrt(n))   [SOL22 PEI term] 

       CALL ReadI(line,lambda_opt,0,3,'Coulomb logarithm calc opt') 

    else if (tag(1:3) == 'T48') then 

       !    T48 Coulomb logarithm calculation options 
       !     Coulomb logarithm constant value - default value is 15.0 - this allows 
       !     specification of alternate constant values for option 0. 

       CALL ReadR(line,lambda_val,0.0,HI,'Coulomb logarithm const val') 

    else if (tag(1:3) == 'T49') then 

       ! T49 Blob frequency for pinch velocity. 
       ! 
       ! T29 is the PDF of the radial velocities. An additional option 
       ! is to specify the rate of blobs for the entire plasma. This is 
       ! to effect that the transport will be intermittent. I.e., instead 
       ! of always sampling the velocity PDF at every step, sample it 
       ! proportional to the amount of blobs it would see. For instance, 
       ! if fblob = 1e6 and qtim = 1e-8, then then probability of choosing 
       ! a velocity from the PDF is: 
       ! prob_choosing = fblob * qtim = 1e-2. 
       ! In words, the ion sees 0.01 blobs every time step. So the 
       ! probability of choosing from the PDF is 1%, which is easily done. 
       ! If fblob = 0.0 this option has no effect. 

       call readr(line, fblob, -1.0, HI, 'Blob frequency') 

    else if (tag(1:3) == 'T50') then 

       ! T50 Core diffusion coefficient. 
       ! -1.0 = Same as CDPERP. 
       call readr(line, cdperpc,-1.0,HI, 'Core diffusion coefficient') 

    else if (tag(1:3) == 'T54') then 

       ! T51-53 are free to use. They were old options I scrapped. 

       call readi(line, balloon_opt, 0, 1, &
            'Ballooning transport approx. switch') 

    else if (tag(1:3) == 'T55') then 

       ! T55  Divertor radial velocity factor for PDF option. If the ion 
       ! is below/above the X point (LSN/USN), then multiply the radial 
       ! velocity when chosen by this value. This is to simulate lower 
       ! (<1) or higher (>1) radial transport in the divertor. Experiments 
       ! have shown lower before. 

       call readr(line, div_vr_fact, -HI, HI, &
            'Divertor radial velocity factor') 

    else if (tag(1:3) == 'T56') then 
       call readi(line, in_blob_switch, 0, 1, &
            'Turn off parallel transport when in blob switch') 
    else if (tag(1:3) == 'T57') then 

       ! T57 T58 
       ! Options for checking inward moving hole-like transport near 
       ! separatrix. 

       call readi(line, hole_switch, 0, 1, &
            'Check for holes near separatrix switch') 
    else if (tag(1:3) == 'T58') then 
       call readr(line, hole_lambda, 0.0, HI, &
            'Hole frequency decay length') 

    else if (tag(1:3) == 'T59') then 
       ! T59 Core pinch value. 
       call readr(line, core_pinch, -HI, HI, &
            'Core pinch value') 

    else if (tag(1:3) == 'T60') then 
       ! T60 Minium psin for blob-like impurity transport model. 
       call readr(line, blob_min_rmrsomp, -HI, HI, &
            'Minimum R-Rsep @ OMP for blob-like transport model') 

    else if (tag(1:3) == 'T61') then 
       ! T61 Birth location of blob/holes. 
       call readr(line, blob_birth_rmrsomp, -HI, HI, &
            'R-Rsep @ OMP for blob/hole birth') 

    else if (tag(1:3) == 'T62') then 
       ! T62 exponential decay length for blob frequency. 
       call readr(line, blob_lambda, -HI, HI, &
            'Blob frequency decay length') 



       
       !===================================================
       !
       ! TAG Series V
       !
       !===================================================
       ! ----------------------------------------------------------------------- 
       !...  Vacuum grid options: 
    else if (tag(1:3) == 'V01') then 
       if (line(7:9) == '1.0') then 
          CALL RDRARN(vacseg,vacnseg,MAXNKS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,3,'Vacuum grid data 1.0',IERR) 

          DO i1 = 1, vacnseg 
             ! jdemod - changed to f12 due to some formatting errors in INTEL compiler 
             WRITE(fp,'(A,1P,4G12.4)') 'VS: ',(vacseg(i1,i2),i2 = 1,4) 
          end do

          !          iflexopt(2) = 20 
       ELSE 
          CALL ER('RUI','Unsupported version for V01 VACSEG',*99) 
       end if

    else if (tag(1:3) == 'V02') then 
       if     (line(7:9) == '1.3') then 
          CALL RDRARN(vacpla,vacnpla,MAXNKS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,7,'Vacuum grid plasma 1.3',IERR) 
          DO i1 = 1, vacnpla 
             if (vacpla(i1,1) /= -3.0.and.vacpla(i1,1) /= -4.0.and. &
                  vacpla(i1,1) /=  0.0 ) &
                  CALL ER('UnstrInput','V02 1.3 must have type -3.0, -4.0 '// &
                  'or 0.0',*99) 
          end do
          DO i1 = 1, vacnpla 
             ! jdemod - changed to g12 due to some formatting errors in INTEL compiler 
             WRITE(fp,'(A,1P,8G12.4)') 'VS: ',(vacpla(i1,i2),i2 = 1,8) 
          end do
       else if (line(7:9) == '1.2') then 
          CALL RDRARN(vacpla,vacnpla,MAXNKS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,6,'Vacuum grid plasma 1.2',IERR) 
          DO i1 = 1, vacnpla 
             ! jdemod - changed to g12 due to some formatting errors in INTEL compiler 
             WRITE(fp,'(A,1P,7g12.4)') 'VS: ',(vacpla(i1,i2),i2 = 1,7) 
          end do
       else if (line(7:9) == '1.1') then 
          CALL RDRARN(vacpla,vacnpla,MAXNKS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,7,'Vacuum grid plasma 1.1',IERR) 
          DO i1 = 1, vacnpla 
             ! jdemod - changed to g12 due to some formatting errors in INTEL compiler 
             WRITE(fp,'(A,1P,8g12.4)') 'VS: ',(vacpla(i1,i2),i2 = 1,8) 
          end do
       else if (line(7:9) == '1.0') then 
          CALL RDRARN(vacpla,vacnpla,MAXNKS,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,6,'Vacuum grid plasma 1.0',IERR) 
          !...      Adjust input to make compatible with version 1.1: 
          DO i1 = 1, vacnpla 
             DO i2 = 8, 3, -1 
                vacpla(i1,i2) = vacpla(i1,i2-1) 
             end do
             vacpla(i1,2) = vacpla(i1,1) 
          end do
          DO i1 = 1, vacnpla 
             ! jdemod - changed to g12 due to some formatting errors in INTEL compiler 
             WRITE(fp,'(A,1P,8g12.4)') 'VS: ',(vacpla(i1,i2),i2 = 1,8) 
          end do
       ELSE 
          CALL ER('RUI','Unsupported revision for V02 VACPLA',*99) 
       end if



       !===================================================
       !
       ! TAG Series W
       !
       !===================================================

    else if (tag(1:3) == 'w01') then 
       !     TAGS W01 and W02 

       !     Options for walls 

       !     wall_plasma_opt - option to select algorithm to calculate 
       !     plasma conditions associated with wall elements 

       CALL ReadI(line,wall_plasma_opt,0,2,'Wall Plasma Option') 
    else if (tag(1:3) == 'w02') then 

       !     wall_plasma_fact - scaling factor used by wall_plasma_opt 
       !     algorithms 

       CALL ReadR(line,wall_plasma_fact,0.0,HI, &
            'Wall Plasma Scaling Factor') 

    ELSE 
       CALL ER('ReadUnstructuredInput','Unrecognized tag',*99) 
    end if



    !write(0,*) 'RUI END  :',trim(tag),':',trim(line),':'

    
    !...  Check tag list: 
    DO i1 = 1, ntaglist-1 
       if (tag == taglist(i1)) then 
          CALL ER('ReadUnstructuredInput','Duplicate tag detected',*10) 
10        WRITE(0,*) 'TAG: "'//tag//'"' 
          WRITE(0,*) 'PROGRAM STOP' 
          STOP 
       end if
    end do

    RETURN 

    ! There is an error: 

99  WRITE(SLOUT,'(5X,3A)') 'LINE = "',line,'"' 
    WRITE(SLOUT,'(5X,3A)') 'TAG  = "',tag ,'"' 
    WRITE(0    ,'(5X,3A)') 'LINE = "',line(1:LEN_TRIM(line)),'"' 
    WRITE(0    ,'(5X,3A)') 'TAG  = "',tag ,'"' 
    WRITE(0,*) '    DIVIMP HALTED' 
    STOP 
  END subroutine readunstructuredinput

  ! ====================================================================== 

  REAL FUNCTION GetInputVersion(line) 
    IMPLICIT none 

    CHARACTER*(*) line 

    integer :: i,j,k 

    GetInputVersion = -1.0 

    !...  Find location of first 2 spaces on the line: 
    i = 0 
    j = 0 
    DO k = 1, LEN_TRIM(line) 
       if (line(k:k) == ' '.and.i /= 0) j = k 
       if (line(k:k) == ' '.and.i == 0) i = k 
       if (j /= 0) exit 
    end do

    if (i /= 0.and.k /= 0.and.j-1 > 1) &
         READ(line(i+1:j-1),*) GetInputVersion 

    RETURN 
    ! 99   STOP 
  END function getinputversion


  SUBROUTINE ReadTagSeries_G(tag,line,fp) 
    use subgrid_options 
    use ribbon_grid_options 
    use mod_params 
    use mod_slcom 
    use mod_comtor 
    use mod_cgeom 
    use mod_io 
    IMPLICIT none 

  ! ====================================================================== 

  ! subroutine: ReadTagSeries_G 

  ! Unstructure input related to the magnetic grid are loaded.  The input 
  ! tag starts with 'G'. 

    !     READ "G" Series Unstructured input 


    integer :: fp 
    character :: line*(*),tag*3 

    integer :: ierr,i1,i2,ir 

    !
    !  Inputs for TAG series G
    !
    if (tag(1:3) .eq. 'G01') then
       !   Sample input 
       !   '+G01    Grid Optgion    0-JET 1-ASDEX 2-ITER             '     3
       call divrd(cgridopt,.true.,0,.true.,9,'GRID OPTION          ',ierr)
    elseif (tag(1:3) .eq. 'G02') then
       !   Sample input 
       !   '+G02    Non-Orthogonal Grid option 0-off 1-JET N.O.     '     3
       call divrd(northopt,.true.,0,.true.,3,'NON-ORTHO. GRID OPT  ',ierr)
    elseif (tag(1:3) .eq. 'G03') then
       !   Sample input 
       !   '+G03    Parallel Distance Option 0-centers 1-edges      '     1
       call divrd(pdopt,.true.,0,.true.,1,   'PARALLEL DIST (S) OPT',ierr)
    elseif (tag(1:3) .eq. 'G04') then
       !   Sample input 
       !   '+G04    Cross-field Distance Option 0-centres 1-edges   '     1
       call divrd(cfdopt,.true.,0,.true.,1,'Proper Cross-field dist',ierr)
    elseif (tag(1:3) .eq. 'G05') then
       !   Sample input 
       !   '+G05    RZ calculation option 0-centers 1-Actual RZ     '     2
       call divrd(rzopt,.true.,0,.true.,3,   'CALCULATE ACTUAL R,Z ',ierr)
    elseif (tag(1:3) .eq. 'G06') then
       !   Sample input 
       !   '+G06    XY Grid option 0-off 1-on                       '     0
       call divrd(xygrid,  .true.,0,.true.,1,'XY GRID OPTION       ',ierr)
    elseif (tag(1:3) .eq. 'G07') then
       !   Sample input 
       !   '+G07    Cell Area Calculation Option 0-approx 1-polygon '     1
       call divrd(cvolopt,.true.,0,.true.,1, 'Cell Volumes from PGs',ierr)
    elseif (tag(1:3) .eq. 'G08') then
       !   Sample input 
       !   '+G08 T   Ion Wall Option        0 to 2                   '     2
       call divrd(CIONR ,.TRUE., 0,.TRUE., 2,'ION WALL OPTION      ',IERR)
    elseif (tag(1:3) .eq. 'G09') then
       !   Sample input 
       !   '+G09 T   Neutral Wall Option    0 to 4                   '     4
       call divrd(CNEUR ,.TRUE., 0,.TRUE., 7,'NEUTRAL WALL OPTION  ',IERR)
    elseif (tag(1:3) .eq. 'G10') then
       !   Sample input 
       !   '+G10 T   Trap Wall Option       0 to 4                   '     4
       call divrd(CTRAP ,.TRUE., 0,.TRUE., 8,'TRAP WALL OPTION     ',IERR)
    elseif (tag(1:3) .eq. 'G11') then
       !   Sample input 
       !   '+G11 T   Vessel Wall Redefinition Option (Baffle Removal)'     0
       call divrd(redefopt,.true.,0,.true.,1,'VESSEL REDEF. OPT    ',IERR)
    elseif (tag(1:3) .eq. 'G12') then
       !   Sample input 
       !   '+G12 T   Target Position Option 0 to 6                   '     6
       call divrd(CTARGOPT,.TRUE.,0,.TRUE.,6,'TARGET POSITION OPT  ',IERR)
    elseif (tag(1:3) .eq. 'G13') then
       !   Sample input 
       !   '+G13 T   Pre-defined geometry option -1-off 0-719 1-307  '    -1
       call divrd(CGEOOPT,.TRUE.,-1,.TRUE.,1,'GEOMETRY OPTION      ',IERR)
    elseif (tag(1:3) .eq. 'G14') then
       !   Sample input 
       !   '+G14    Central Mirror Ring Location       (IR)         '     1
       call divrd(ircore,.true.,1,.false.,0,'CORE MIRROR RING SPEC' ,ierr)
    elseif (tag(1:3) .eq. 'G15') then
       !   Sample input 
       !   '+G15    Rectangular grid for neutrals 0calculate 99file '      0
       call divrd(CRECT ,.TRUE. , 0 ,.TRUE. , 99,'RECT GRID CALC.   ',IERR)
    elseif (tag(1:3) .eq. 'G16') then
       !   Sample input 
       !   '+G16 ' 'TN    Set of Plate coordinates                  '
       !   '    TN    Ring #, Outer R,Z   , Inner R,Z :      '       0
       call divrd(PLATCO,NPLAT,MAXNRS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,4,'PLATE COORDINATES',IERR)
    elseif (tag(1:3) .eq. 'G17') then
       !   Sample input 
       !   '+G17 ' 'Wall coordinates                                '
       !   '    TN    R,Z coordinates starting at outer plate   '   0
       call divrd(WALLCO,NWALL,MAXPTS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,1,'WALL COORDINATES',IERR)
    elseif (tag(1:3) .eq. 'G18') then
       !   Sample input 
       !   '+G18 ' 'Wall coordinates - PFZ                          '
       !   '    TN    R,Z coordinates                           '  0
       call divrd(WALLCO2,NWALL2,MAXPTS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,1,'WALL COORDINATES',IERR)
    elseif (tag(1:3) .eq. 'G19') then
       !   Sample input 
       !   '+G19 ASDEX U - GRID CHARACTERISTICS:  Number of Rings    '    26
       !
       !     NOTE: These items are now usually on the GEOM: line at the beginning
       !           of the grid file or are calculated automatically from the grid
       !     The following set of numbers specify the characteristic values of
       !     the ASDEX U style grid file to be read in. The Number of rings
       !     and elements per ring (a ring is a set of polygons in a row or
       !     indexed by the same ring number), the cutring,(which specifies the
       !     end of  rings for the core and trapped plasma), and the cutpts 1
       !     and 2 which specify the points on the rings numbered less than the
       !     cut ring where the splits for core and trapped or private plasma
       !     occur.
       !
       call divrd(maxrings,.TRUE.,0  ,.true.,maxnrs,'MAXRINGS in AUG',IERR)
    elseif (tag(1:3) .eq. 'G20') then
       !   Sample input 
       !   '+G20                                 Number of Knots    '    34
       call divrd(maxkpts,.TRUE.,0  ,.true.,maxnks,'MAX PTS in AUG',IERR)
    elseif (tag(1:3) .eq. 'G21') then
       !   Sample input 
       !   '+G21                                 Cut ring           '     7
       call divrd(cutring,.TRUE.,0  ,.true.,maxrings,'CUTRING in AUG',IERR)
    elseif (tag(1:3) .eq. 'G22') then
       !   Sample input 
       !   '+G22                                 Cut point 1        '     1
       call divrd(cutpt1,.TRUE.,0  ,.true.,maxkpts+1,'CUTPT1 in AUG',IERR)
    elseif (tag(1:3) .eq. 'G33') then
       !   Sample input 
       !   '+G33                                 Cut point 2        '    34
       call divrd(cutpt2,.TRUE.,0  ,.true.,maxkpts+1 ,'CUTPT2 in AUG',IERR)


       !
       ! jdemod - this option does not appear to be in use and tag collides with existing
       ! 
       !       elseif     (tag(1:3) == 'G01') then 
       !          !...    Not sure if Dave has reserved G01 already, and the web server is down: 
       !          CALL ReadR(line,grd_thresh,0.0,100.0,'Refinement threshold') 

    else if (tag(1:3) == 'G23') then 
       ! ----------------------------------------------------------------------- 

       !     TAG G23: SONNET grid SUB type option 

       !     This option is used to tag a SONNET style grid as being an 
       !     FRC (Field Reversed Configuration) custom grid. At the moment 
       !     only one option is supported but further subtypes could be 
       !     added for Sonnet grids requiring special processing. 

       !     0 = stanard Sonnet grid = Default 
       !     1 = FRC version 1 - type of Sonnet grid 
       !         - used to set various FRC related options 
       !     2 = sonnet grid without boundary cells - boundary cells are added 
       !         - useful for carre grids 

       CALL ReadI(line,sonnet_grid_sub_type,0,2,'SONNET Grid SUB-Type Option') 

    else if (tag(1:3) == 'G34') then 
       ! ----------------------------------------------------------------------- 

       !     TAG G34 : Machine Type 

       !       0 = jet 
       !       1 = diiid 
       !       2 = alcator (cmod) 
       !       3 = aug (asdex upgrade) 
       !       4 = iter 
       !       5 = ignitor 
       !       6 = fire 


       CALL ReadI(line,tmachine_opt,0,6,'MACHINE TYPE option') 

    else if (tag(1:3) == 'G35') then 

       ! ----------------------------------------------------------------------- 

       !     TAG G35 : Shot Identification string for JET cataloguing system 


       !       Shot identification string for cataloguing on JET 

       CALL ReadC(line,divshotid,'OPTIONAL SHOT ID') 

    else if (tag(1:3) == 'G36') then 
       ! ----------------------------------------------------------------------- 

       !     TAG G36 : Parallel Ion Reflection option 


       !       Option to activate parallel ion reflection 

       CALL ReadI(line,s_reflect_opt,0,1,'PARALLEL ION REFLECTION OPT') 

    else if (tag(1:3) == 'G37') then 
       ! ----------------------------------------------------------------------- 

       !...    Data for modifying (cutting/extending) the magnetic grid: 

       if (line(7:9) == '1.0') then 
          !...      Inputs: 

          !         Type - 1: Cut ring to specified line segment. 
          !                2: Extend end of ring to specified line segment. 
          !         Mode - 1: Work from low index target (and remove portion from cut to outer target if 
          !                   necessary). 
          !                2: Work from high index target (and remove portion from cut to outer target if 
          !                   necessary) 
          !         Refinement - 0: none 
          !         Range1 - Start of IR index range of rings to be cut (inclusive). 
          !         Range2 - End of IR index range of rings to be cut. 

          CALL RdRarn(grdmod,grdnmod,1000,-MACHHI,MACHHI,.FALSE., &
               -MACHHI,MACHHI,4,'Magnetic grid modification',ierr) 
          if (ierr /= 0) call er('readtagseriesg','error grdmod',*99) 

          WRITE(fp,*) 
          WRITE(fp,'(A)') 'Grid modification data:' 
          DO i1 = 1, grdnmod 
             WRITE(fp,'(5X,5F9.5)') (grdmod(i1,i2),i2 = 1,5) 
          end do

       end if

    else if (tag(1:3) == 'G38') then 
       ! ----------------------------------------------------------------------- 

       !     The following tags are related to the subgrid option for recording 
       !     more detailed data on a finer grid 

       !     G38: Base subgrid option ON/OFF 
       !     G39: R,Z dimensions of the grid region 
       !     G40: RMIN,RMAX of the subgrid region 
       !     G41: ZMIN,ZMAX of the subgrid region 

       ! ----------------------------------------------------------------------- 

       !     TAG G38 : Subgrid option 


       !     G38: Option to activate the subgrid data collection 

       CALL ReadI(line,subgrid_opt,0,1,'BASE SUBGRID OPTION') 

    else if (tag(1:3) == 'G39') then 

       !     G39: R,Z dimensions of the grid region 

       CALL Read2I(line,sg_rdim,sg_zdim,1,500,'SUBGRID RDIM,ZDIM') 

    else if (tag(1:3) == 'G40') then 

       !     G40: RMIN,RMAX of the subgrid region 

       CALL Read2R(line,sg_rmin,sg_rmax,-HI,HI,'SUBGRID RMIN,RMAX') 

    else if (tag(1:3) == 'G41') then 

       !     G41: ZMIN,ZMAX of the subgrid region 

       CALL Read2R(line,sg_zmin,sg_zmax,-HI,HI,'SUBGRID ZMIN,ZMAX') 

    else if (tag(1:3) == 'G42') then 

       ! ----------------------------------------------------------------------- 

       !     Options related to ribbon grids 
       !     G42 - grid generation option - <i4> 
       !     G43 - intersection point averaging option - opt_block_av - <r4> 
       !     G44 - maximum R separation in grid generator - max_r_sep - <r4> 
       !     G45 - maximum S/Z separation in grid generator - max_s_sep - <r4> 
       !     G46 - min number of cells on ring - min_cells - <i4> 
       !     G47 - castem output identifier - <string> 
       !     G48 - min and max S for selecting intersection subset  2 x <r4> 
       !     G49 - min and max R for intersection subset grid generation 2 x <r4> 
       !     G50 - min and max S for intersection subset grid generation 2 x <r4> 

       !     G42: Ribbon grid option  0=unstructured  1=structured 

       CALL ReadI(line,rg_grid_opt,0,1,'RIBBON GRID OPTION') 
    else if (tag(1:3) == 'G43') then 

       !     G43: Ribbon grid option  1=unstructured  2=structured 

       CALL ReadI(line,rg_block_av,0,1,'BLOCK AVERAGE OPTION') 

    else if (tag(1:3) == 'G44') then 

       !     G44: Maximum row separation 

       Call ReadR(line,rg_max_r_sep,0.0,HI,'Maximum Row separation (m)') 
    else if (tag(1:3) == 'G45') then 

       !     G45: Maximum cell separation 

       Call ReadR(line,rg_max_s_sep,0.0,HI,'Maximum cell separation (m)') 

    else if (tag(1:3) == 'G46') then 

       !     G46: Minimum number of cells in a row 

       CALL ReadI(line,rg_min_cells,1,9999,'Minimum cells in a row') 

    else if (tag(1:3) == 'G47') then 

       !     G47: Castem data set identifier 

       CALL ReadC(line,rg_castem_data,'CASTEM DATA SET IDENTIFIER') 
       write(0,*) 'RG_CASTEM_DATA:',trim(rg_castem_data),':' 

    else if (tag(1:3) == 'G48') then 

       !     G48 - min and max S for selecting intersection subset  2 x <r4> 

       CALL Read2R(line,rg_int_win_mins,rg_int_win_maxs,-HI,HI,'RIBBON GRID INTERSECTION SUBSET RANGE [S1,S2]') 

    else if (tag(1:3) == 'G49') then 

       !     G49 - min and max S for selecting intersection subset  2 x <r4> 

       CALL Read2R(line,rg_minr,rg_maxr,-HI,HI,'RIBBON GRID SUBSET R RANGE [R1,R2]') 

    else if (tag(1:3) == 'G50') then 

       !     G50 - min and max S for intersection subset grid generation 2 x <r4> 

       CALL Read2R(line,rg_mins,rg_maxs,-HI,HI,'RIBBON GRID SUBSET S RANGE [S1,S2]') 

    else if (tag(1:3) == 'G51') then 

       !     G51 - cutoff factor for ring generation - rings 
       !           in a ribbon grid with a length factor smaller than 
       !           this value will not be generated. 

       CALL ReadR(line,lcutoff,-HI,HI,'RING CUTOFF LENGTH FACTOR') 

    else if (tag(1:3) == 'G52') then 

       !     G52 - Cell spacing option ... option to calculate the cell 
       !           boundary spacing along a ring. Only option 0 is currently 
       !           available which uses an exponential factor given in G53. 
       !           A cell_spacing_factor of 1 gives a linear spacing 

       CALL ReadI(line,cell_spacing_option,0,0,'CELL SPACING OPTION') 

    else if (tag(1:3) == 'G53') then 

       !     G53 - cell spacing factor 
       !           Used to determine cell boundary spacing along the rings 

       CALL ReadR(line,cell_spacing_factor,-HI,HI,'CELL SPACING FACTOR') 

    else if (tag(1:3) == 'G54') then 

       !     G54 - Cell spacing option ... option to calculate the cell 
       !           boundary spacing along a ring. Only option 0 is currently 
       !           available which uses an exponential factor given in G53. 
       !           A cell_spacing_factor of 1 gives a linear spacing 

       CALL ReadI(line,ribbon_input_format_opt,0,1, &
            'INTERSECTION DATA INPUT DATA FORMAT OPTION: 0=CASTEM 1 = RAY') 


    else if (tag(1:3) == 'G55') then 

       ! ----------------------------------------------------------------------- 
       !     G55 - read in a list of ikoffsets to move the center of the 
       !           ring for background plasma calculation 

       !           n_ikoffsets 

       !           option   ir1   ir2    offset_value 

       !           option=0 offset = ik_mids + offset_value (index offset) 
       !           option=1 offset = offset_value * SMAX    (fractional offset) 
       !           option=2 offset = offset_value * PMAX    (fractional offset) 


       CALL RDRARN(ik_offset_data,n_ik_offsets,MAXNRS, &
            -MACHHI,MACHHI,.FALSE., &
            -MACHHI,MACHHI,3,'IK OFFSET DATA',IERR) 


    else if (tag(1:3) == 'G56') then 

       !     jdemod - sol22_halfringlen_opt = 0 = ringlen/2 
       !                                    = 1 = ksb(ikmid,ir) 
       !            - the ikmid value should be at ~ 1/2 the ring length 
       !            - note: using option 1 will support calculations that 
       !                    move the midpoint 

       CALL ReadI(line,sol22_halfringlen_opt,0,1,'SOL22 half ring length determination option') 

    ELSE 
       CALL ER('ReadTagSeriesG','Unrecognized tag',*99) 
    end if

    RETURN 
    !.... Error: 
99  WRITE(SLOUT,'(5X,3A)') 'LINE = "',line,'"' 
    WRITE(SLOUT,'(5X,3A)') 'TAG  = "',tag ,'"' 
    WRITE(0    ,'(5X,3A)') 'LINE = "',line(1:LEN_TRIM(line)),'"' 
    WRITE(0    ,'(5X,3A)') 'TAG  = "',tag ,'"' 
    WRITE(0,*) '    DIVIMP HALTED' 
    STOP 
  END subroutine readtagseries_g


  SUBROUTINE ReadTagSeries_H(line,tag,fp) 
    Use comhc            ! Assign values to Hydrocarbon following common block. 
    use hc_kinetics_options 
    use mod_comtor 
    use mod_slcom 
    use mod_hc_global_opts 
    use mod_io 
    use mod_pindata
    IMPLICIT none 

  ! ====================================================================== 

  ! subroutine: ReadTagSeries_H 

  ! Unstructured input data tagged with an H are loaded. The majority of 
  ! these refer to hydrocarbon code modeling options. 

    !     READ "H" Series Unstructured input 

    integer :: ierr,in
    integer :: fp 
    character :: line*(*),tag*3 

    !
    !  Inputs for TAG series H
    !
    if (tag(1:3) .eq. 'H01') then
       !   Sample input 
       !   '+H01    PIN Random number seed  (<0=1, 0 generate new)  '      0
       call divrd(PINISEED,.FALSE.,0,.FALSE.,0,'PIN RANDOM NUM. SEED',IERR)
    elseif (tag(1:3) .eq. 'H02') then
       !   Sample input 
       !   '+H02    PIN Data Print option  (0 reduced, 1 more)      '      0
       call divrd(PINPRINT,.true., 0 ,.true.,1,'PIN DATA PRINT OPT',IERR)
    elseif (tag(1:3) .eq. 'H03') then
       !   Sample input 
       !   '+H03 TN408 Run PIN from inside DIVIMP  0-NO 1-YES        '    1
       call divrd(CPINOPT,.TRUE.,0  ,.TRUE. ,4  ,'RUNPIN 0-NO 1-YES',IERR)
    elseif (tag(1:3) .eq. 'H04') then
       !   Sample input 
       !   '+H04 ' 'TN408 Pin: reire07                              '
       !   READ(CPINCOM(11:80),'(A69)') ACTPIN
       call divrd(CPINCOM,'COMMAND TO RUN PIN',IERR)
       ! Use : as delimiter if present - otherwise first 10 characters go to line label. 
       in = index(cpincom,':')
       if (in.gt.0) then 
          READ(CPINCOM(in+1:),'(a69)') ACTPIN
       else
          READ(CPINCOM(11:),'(a69)') ACTPIN
       endif
       !write(0,*) 'cpincom:',trim(cpincom),':',in,':',trim(cpincom(in+1:)),':',trim(actpin),':'
    elseif (tag(1:3) .eq. 'H05') then
       !   Sample input 
       !   '+H05      PIN Cell Area Option (IHCORR)                 '    1
       call divrd(IHCORR,.TRUE.,0  ,.TRUE.,1 ,'PIN CELL AREA OPTION',IERR)
    elseif (tag(1:3) .eq. 'H06') then
       !   Sample input 
       !   '+H06      PIN Hybrid Wall Option 0=off 1,2=selection    '    0
       call divrd(IHYBRID,.TRUE.,0  ,.TRUE. ,6 ,'HYBRID WALL IN PIN',IERR)
    elseif (tag(1:3) .eq. 'H07') then
       !   Sample input 
       !   '+H07      PIN Puffing Option - 0=off 1=on               '    0
       call divrd(pinpuff,.true.,0,.true.,2,    'PIN PUFFING OPTION',ierr)
    elseif (tag(1:3) .eq. 'H08') then
       !   Sample input 
       !   '+H08      PIN Puff Location switch - 0=main SOL 1=PP    '    0
       call divrd(swpvhpf,.true.,0,.true.,1,  'PUFF LOCATION OPTION',ierr)
    elseif (tag(1:3) .eq. 'H09') then
       !   Sample input 
       !   '+H09      PIN Puff fraction (opt 1)                     '   0.0
       call divrd(hpcpuf, .TRUE., 0.0,.FALSE.,0.0,'PIN RE-PUFF FRAC',IERR)
    elseif (tag(1:3) .eq. 'H10') then
       !   Sample input 
       !   '+H10      PIN Recycle -> Puff fraction (puff opt 2)     '   0.16
       call divrd(ppcpuf,.true.,0.0,.FALSE.,0.0, 'FLUX-> PUFF OPT 2',IERR)
    elseif (tag(1:3) .eq. 'H11') then
       !   Sample input 
       !   '+H11      PIN Puff Injection temperature (eV)           '   0.5
       call divrd(tpufh, .TRUE., 0.0,.FALSE.,0.0,'PIN RE-PUFF TEMP ',IERR)
    elseif (tag(1:3) .eq. 'H12') then
       !   Sample input 
       !   '+H12      PIN Puff location indices JHPUF1(1 and 2)'  -17 -1000
       call divrd(jhpuf1(1),jhpuf1(2),.false.,0,.false.,1,'PUFF LOCATION INDICES 1',ierr)
    elseif (tag(1:3) .eq. 'H13') then
       !   Sample input 
       !   '+H13      PIN Puff location indices JHPUF2(1 and 2)'  -16 -1001
       call divrd(jhpuf2(1),jhpuf2(2),.false.,0,.false.,1,'PUFF LOCATION INDICES 2',ierr)


    elseif (tag(1:3) == 'H15') then 
       ! ammod begin. 
       ! ---------------------------- 

       ! Options added for hydrocarbon following. 

       ! ----------------------------

       Call ReadI(line,hc_follow_option,0,1,'Hydrocarbon following option, 0-off, 1-on') 

       !       Set global_hc_follow_option to match hc value read in 
       !       Initialization of global_hc_follow_option occurs in the setup 
       !       routine where initialization of the hc_follow_option also takes place 

       global_hc_follow_option = hc_follow_option 

    else if (tag(1:3) == 'H16') then 
       Call ReadI(line,hc_higher_hcs_option,0,1,'Follow higher hydrocarbon (C2+) option, 0-off, 1-on') 
    else if (tag(1:3) == 'H17') then 
       Call ReadI(line,hc_wbc_comp_option,0,1,'WBC comparison case, 0-off, 1-on') 

    else if (tag(1:3) == 'H20') then 
       Call ReadI(line,hc_sputtering_model,0,1, &
            'Model for sputtered species release, 0-preset, 1-Mech') 
    else if (tag(1:3) == 'H21') then 
       Call ReadI(line,hc_sputtered_hc_species,0,58,'Preset sputtered hydrocarbon species, 10-Methane(CH4)') 
    else if (tag(1:3) == 'H22') then 
       Call ReadI(line,hc_evolution_model_primary,1,3,'Model for HC data primary, 1-E&L, 2-Brooks, 3-Janev') 
    else if (tag(1:3) == 'H23') then 
       Call ReadI(line,hc_evolution_model_secondary,0,3,'Model for HC data secondary, 0-none, 1-E&L, 2-Brooks,3-Janev') 
       !     Available libraries 1) Ehrhardt and Langer (PPPL, 1987) 
       !                         2) Alman, Ruzic, Brooks (Phys. Plasmas, 2000) 
       !                         3) Janev, et al. (NIFS, 2001) 
    else if (tag(1:3) == 'H24') then 
       Call ReadI(line,hc_launch_location,-1,6,'Model for HC launch location (same options as CNEUTB)') 

       !        Values assigned in global_hc_assign_inputs after the 
       !        entire input file has been read in 

       !	If (hc_launch_location .eq. -1) Then 
       !		hc_launch_location = CNEUTB 
       !	End If 
    else if (tag(1:3) == 'H25') then 
       Call ReadI(line,hc_launch_angle_velocity,-1,20,'Model for HC launch angle/velocity (same options as CNEUTC)') 

       !        Values assigned in global_hc_assign_inputs after the 
       !        entire input file has been read in 

       !	If (hc_launch_angle_velocity .eq. -1) Then 
       !		hc_launch_angle_velocity = CNEUTC 
       !	End If 
    else if (tag(1:3) == 'H26') then 
       Call ReadI(line,hc_launch_velocity_model,0,2,'Launch velocity model, 0-const, 1-MB dist., 2-dual MB') 
    else if (tag(1:3) == 'H27') then 
       Call ReadR(line,hc_dual_mb_pri_vel_flux,0.0,1.0,'Dual MB velocity flux primary MB contrib., 0.0-1.0') 
    else if (tag(1:3) == 'H28') then 
       Call ReadR(line,hc_dual_mb_sec_mean_temp,0.0,2000.0,'Dual MB launch velocity T2 (deg K), 0.0-2000.0') 

    else if (tag(1:3) == 'H30') then 
       Call ReadI(line,hc_neut_ion_velocity,-1,3,'Neutral->Ion initial velocity (same options as CNEUTG)') 

       !        Values assigned in global_hc_assign_inputs after the 
       !        entire input file has been read in 

       !       If (hc_neut_ion_velocity .eq. -1) Then 
       !		hc_neut_ion_velocity = CNEUTG 
       !	End If 
    else if (tag(1:3) == 'H31') then 
       Call ReadI(line,hc_ion_neut_angle,0,2,'Ion->neutral angle emission, 0-isotropic, 1-sine, 2-S dir') 
    else if (tag(1:3) == 'H32') then 
       Call ReadI(line,hc_ion_neut_velocity,0,1,'Ion->neutral velocity, 0-ion energy, 1-') 
    else if (tag(1:3) == 'H33') then 
       Call ReadI(line,hc_lambda_calc,0,1,'Improved calculation for lambda (Sivukhin),0-off,1-on') 
    else if (tag(1:3) == 'H34') then 
       Call ReadI(line,hc_disable_transitions,0,1,'Disable HC transitions, 0-off,1-on') 
    else if (tag(1:3) == 'H35') then 
       Call ReadI(line,hc_presheath_efield,0,1,'Improved model for electric field force,0-off,1-on') 
    else if (tag(1:3) == 'H36') then 
       Call ReadR(line,hc_efield_drop_fraction,0.0,1.0,'Fraction of potential drop in Debye region,0.0-1.0') 
    else if (tag(1:3) == 'H37') then 
       Call ReadI(line,hc_efield_cells,0,5,'Cells from target to apply improved e-field, 0-5') 

    else if (tag(1:3) == 'H40') then 
       Call ReadI(line,hc_neutral_reflection_option,0,1,'Neutral HC reflection switch') 
    else if (tag(1:3) == 'H41') then 
       Call ReadI(line,hc_ion_reflection_option,0,1,'Ion HC reflection switch') 
    else if (tag(1:3) == 'H42') then 
       Call ReadI(line,hc_reflection_coef_model,0,4,'Reflection model, 0-preset, 1-Janev, 2-A&R, 3-CH4 = 1.0') 
    else if (tag(1:3) == 'H43') then 
       Call ReadR(line,hc_reflection_coef_preset,0.0,1.0,'Preset reflection coef') 
    else if (tag(1:3) == 'H44') then 
       Call ReadI(line,hc_reflection_species_model,0,1,'Reflected species model, 0-preset, 1-Alman&Ruzic') 
    else if (tag(1:3) == 'H45') then 
       Call ReadI(line,hc_reflection_energy_model,0,3,'Reflection energy model,0-set,1-impact,2-thermal,3-AR') 
    else if (tag(1:3) == 'H46') then 
       Call ReadR(line,hc_refl_energy_neutral_preset,0.0,1E4,'Preset reflecting particle energy, neutral impact (eV)') 
    else if (tag(1:3) == 'H47') then 
       Call ReadR(line,hc_refl_energy_ion_preset,0.0,1E4,'Preset reflecting particle energy, ion impact (eV)') 
    else if (tag(1:3) == 'H48') then 
       Call ReadI(line,hc_reflection_angle_model,-1,11,'Refl. angle model,-1 = NRFOPT,0-off,1-spec,2-isotr,3-norm,4-A&R ') 

       !        Values assigned in global_hc_assign_inputs after the 
       !        entire input file has been read in 

       !        If (hc_reflection_angle_model .eq. -1) Then 
       !		hc_reflection_angle_model = NRFOPT 
       !        End If 

    else if (tag(1:3) == 'H50') then 
       Call ReadI(line,hc_sputtering_option,0,1,'HC sputtering switch') 
    else if (tag(1:3) == 'H51') then 
       Call ReadI(line,hc_sticking_coef_model,0,2,'Sticking model, 0-preset, 1-Janev, 2-Alman&Ruzic') 
    else if (tag(1:3) == 'H52') then 
       Call ReadR(line,hc_sticking_coef_preset,-1.0,1.0,'Preset sticking coefficient, -1.0 = CTRESH, 0.0-1.0 stuck') 
    else if (tag(1:3) == 'H53') then 
       Call ReadI(line,hc_sputtering_species_model,0,2,'Sputtered species model, 0-preset, 1-A&R, 2-same') 
    else if (tag(1:3) == 'H54') then 
       Call ReadI(line,hc_sputtering_energy_model,0,3,'Sputtering energy model,0-set,1-impact,2-thermal,3-AR') 
    else if (tag(1:3) == 'H55') then 
       Call ReadR(line,hc_sput_energy_neutral_preset,0.0,1E4,'Preset sputtering particle energy, neutral impact (eV)') 
    else if (tag(1:3) == 'H56') then 
       Call ReadR(line,hc_sput_energy_ion_preset,0.0,1E4,'Preset sputtering particle energy, ion impact (eV)') 
    else if (tag(1:3) == 'H57') then 
       Call ReadI(line,hc_sputtering_angle_model,-1,11, &
            'Sput. angle model,-1 = NRFOPT,0-off,1-spec,2-isosin,' &
            // '3-isocos,4-sqrtcos,5-proj-sqrtsin,10-norm,11-A&R ') 


    else if (tag(1:3) == 'H60') then 
       !--------------------------- 

       ! jdemod 
       !     Wall segment index for HC_launch_location option 6 

       Call ReadI(line,hc_launch6_wall_index,1,maxpts,'Wall segment index for HC launch option 6') 

       !       Set flag indicating that the value has been read 

       hc_launch6_wall_index_set = .true. 

       !       Set flag indicating that the value has been read 

       hc_launch6_wall_index_set = .true. 

    else if (tag(1:3) == 'H61') then 
       !--------------------------- 

       !     HC reaction kinetics option: 
       !     0 = off 
       !     1 = on 
       !     2+= advanced options (to be implemented) 

       Call ReadI(line,hc_reaction_kinetics,0,4,'HC reaction kinetics option') 
    else if (tag(1:3) == 'H62') then 
       !--------------------------- 
       ! jdemod 

       !     HC reaction kinetics option: 
       !     0 = original code 
       !     1 = new implementation of reaction kinetics 

       Call ReadI(line,hc_kinetics_opt,0,1,'NEW HC reaction kinetics option') 


    else if (tag(1:3) == 'H63') then 
       !--------------------------- 
       ! jdemod 

       !     HC kinetic suboption - Vperp evolution: 
       !     0 = Vperp remains constant between reactions 
       !     1 = Vperp = Vpara assigned at reaction start 
       !     2 = Vperp = f(Tperp) - Vperp calculated from temperature 
       !     3 = Vperp diffuses independent of Vpara at same rate 
       Call ReadI(line,hc_vperp_opt,0,3,'HC Perpendicular velocity option') 


    else if (tag(1:3) == 'H64') then 

       !     jdemod - added an input parameter to define the H isotope 
       !              in the hydrocarbon molecules so it can be different 
       !              from the background plasma 

       ! slmod begin - *** TEMP *** 
       Call ReadR(line,input_HC_H_mass,1.0,3.0,'Mass of the H isotope in HC molecules 1.0->3.0') 
       !        STOP 'TURNING OFF FOR NOW...' 
       ! slmod end 


       ! jdemod 
       !--------------------------- 

       !        Values assigned in global_hc_assign_inputs after the 
       !        entire input file has been read in 

       !        If (hc_sputtering_angle_model .eq. -1) Then 
       !		hc_sputtering_angle_model = NRFOPT 
       !	End If 

       ! Hydrocarbon output options. 

    else if (tag(1:3) == 'H90') then 
       Call ReadI(line,hc_coord_print_option,0,1,'Print r,z position data at each timestep, 0-off, 1-on') 
    else if (tag(1:3) == 'H91') then 
       Call ReadI(line,hc_evolve_print_option,0,1,'Print r,z position data at each transition, 0-off, 1-on') 


       ! ---------------------------- 

       ! End hydrocarbon following options. 

       ! ---------------------------- 
       ! ammod end. 

    ELSE 
       CALL ER('ReadTagSeriesH','Unrecognized tag',*99) 
    end if


    RETURN 


99  WRITE(SLOUT,'(5X,3A)') 'LINE = "',line,'"' 
    WRITE(SLOUT,'(5X,3A)') 'TAG  = "',tag ,'"' 
    WRITE(0    ,'(5X,3A)') 'LINE = "',line(1:LEN_TRIM(line)),'"' 
    WRITE(0    ,'(5X,3A)') 'TAG  = "',tag ,'"' 
    WRITE(0,*) '    DIVIMP HALTED' 
    STOP 
  END subroutine readtagseries_h


  SUBROUTINE ReadTagSeries_I(line,tag,fp) 
    use mod_params 
    use mod_comtor 
    use mod_slcom 
    use mod_fperiph_com 
    use mod_cedge2d 
    use mod_promptdep 
    use mod_io 
    use mod_rundiv_local
    IMPLICIT none 

  ! ====================================================================== 

  ! subroutine: ReadTagSeries_I 

  ! Unstructure input related to the magnetic grid are loaded.  The input 
  ! tag starts with 'I'. 

    !     READ "I" Series Unstructured input 

    integer :: ierr
    integer :: fp 
    character :: line*(*),tag*3 

    !
    !  Inputs for TAG series I
    !
    if (tag(1:3) .eq. 'I01') then
       !   Sample input 
       !   '+I01    Injection opt  1/2/3                            '     1
       call divrd(CIOPTE,.TRUE., 0,.TRUE.,14,'INJECTION OPT        ',IERR)
    elseif (tag(1:3) .eq. 'I02') then
       !   Sample input 
       !   '+I02    First diffuse  0inst 1random 2tpara             '     0
       call divrd(CDIFOP,.TRUE., 0,.TRUE., 2,'FIRST DIFFUSE OPT    ',IERR)
    elseif (tag(1:3) .eq. 'I03') then
       !   Sample input 
       !   '+I03    Control switch 0atoms 1ions                     '     0
       call divrd(CNEUTA,.TRUE., 0,.TRUE., 1,'CONTROL SWITCH       ',IERR)
    elseif (tag(1:3) .eq. 'I04') then
       !   Sample input 
       !   '+I04    Self- Sputtering Option 0-off 1-on              '     1
       call divrd(CSELFS,.TRUE., 0,.TRUE., 2,'SELF-SPUTTER OPTION  ',IERR)
    elseif (tag(1:3) .eq. 'I05') then
       !   Sample input 
       !   '+I05    Init ion Vel   1                                '     1
       call divrd(CNEUTG,.TRUE., 0,.TRUE., 3,'INITIAL ION VELOCITY ',IERR)
    elseif (tag(1:3) .eq. 'I06') then
       !   Sample input 
       !   '+I06 TN1465 Follow Imp.Ions Recombined to Neutrals 0=off '     1
       call divrd(CFOLREC,.TRUE., 0,.TRUE.,1,'FOLLOW REC. NEUTRAL  ',IERR)
    elseif (tag(1:3) .eq. 'I07') then
       !   Sample input 
       !   '+I07 TN1479 Ion Prompt Redeposition Option 0=off 1=on    '     0
       call divrd(prompt_depopt,.true.,0,.true.,1,'ION PROMPT DEPOSITION',ierr)
    elseif (tag(1:3) .eq. 'I08') then
       !   Sample input 
       !   '+I08 T   Target Mirror Option 0-off  1-on                '     0
       call divrd(cmiropt, .true.,0,.true.,4,'TARGET MIRROR OPT    ',IERR)
    elseif (tag(1:3) .eq. 'I09') then
       !   Sample input 
       !   '+I09 T   Ion Periphery Option   0 to 3                   '     0
       call divrd(FPOPT,   .TRUE.,0,.TRUE.,6,'ION PERIPHERY OPT    ',IERR)
    elseif (tag(1:3) .eq. 'I10') then
       !   Sample input 
       !   '+I10 TN996 Periphery Recycle Option       0-off 1-on     '     0
       call divrd(fpropt,.true.,  0,.true.,1,'FP RECYCLE OPT       ',ierr)
    elseif (tag(1:3) .eq. 'I11') then
       !   Sample input 
       !   '+I11    Z effective (self)                 Zeff         '      1
       call divrd(CIZEFF,.TRUE. , 0 ,.FALSE., 0 ,'ZEFF(SELF)     ',   IERR)
    elseif (tag(1:3) .eq. 'I12') then
       !   Sample input 
       !   '+I12    Initial ionization state of impurity ions       '      1
       call divrd(CIZSC, .TRUE.,  1, .TRUE.,CION,'INITIAL IZ STATE',  IERR)
    elseif (tag(1:3) .eq. 'I13') then
       !   Sample input 
       !   '+I13    Collision Enhancement Factor       Zenh         '    1.0
       call divrd(CZENH, .FALSE.,0.0,.FALSE.,0.0,'COLLIS ENHANC ZENH',IERR)
    elseif (tag(1:3) .eq. 'I14') then
       !   Sample input 
       !   '+I14    Set Ti = max(Ti,Tb) when reaching state (0 off) '      0
       call divrd(CIZSET,.FALSE., 0 ,.FALSE., 0 ,'SET TI=TB AT STATE',IERR)
    elseif (tag(1:3) .eq. 'I15') then
       !   Sample input 
       !   '+I15    Maximum ionization state                        '      6
       call divrd(NIZS,  .TRUE.,  0, .FALSE.,MAXIZS,'MAX IZ STATE',   IERR)
    elseif (tag(1:3) .eq. 'I16') then
       !   Sample input 
       !   '+I16    Stop following ions reaching Main Plasm 0no 1yes'      0
       call divrd(CSTOP ,.TRUE. , 0 ,.TRUE. , 1 ,'STOP WHEN HIT MAIN',IERR)
    elseif (tag(1:3) .eq. 'I17') then
       !   Sample input 
       !   '+I17    Ion removal loss time              Tloss  (s)   '    0.000
       call divrd(TLOSS,.FALSE.,0.0,.FALSE.,0.0,'ION LOSS TIME', IERR)
    elseif (tag(1:3) .eq. 'I18') then
       !   Sample input 
       !   '+I18 TN480 Ring for ion injection option 2        INJIR  '    1
       call divrd(INJIR,.TRUE.,1,.TRUE.,MAXNRS,'RING NUMBER FOR INJ',IERR)
    elseif (tag(1:3) .eq. 'I19') then
       !   Sample input 
       !   '+I19 TN480 Factor for Region Lower Bound          INJF1  '   0.0
       call divrd(INJF1,.TRUE.,0.0,.TRUE.,1.0  ,'INJECTION AREA LB' ,IERR)
    elseif (tag(1:3) .eq. 'I20') then
       !   Sample input 
       !   '+I20 TN480 Factor for Region Upper Bound          INJF2  '   1.0
       call divrd(INJF2,.TRUE.,INJF1,.TRUE.,1.0,'INJECTION AREA UB' ,IERR)
    elseif (tag(1:3) .eq. 'I21') then
       !   Sample input 
       !   '+I21 TN443 X-max for Far-Periphery region (O/I) '       0.1  0.1
       call divrd(FPXMAXO,fpxmaxI,.TRUE.,0.0,.FALSE.,0.0,'FP DIFFUSION SIZE' ,IERR)
    elseif (tag(1:3) .eq. 'I22') then
       !   Sample input 
       !   '+I22 TN443 Far-Periphery Target Loss Time (O/I) '    1.0e-3  1.0e-3
       call divrd(FPTIMO,fptimI,.TRUE.,0.0,.FALSE.,0.0,'FP TARGET LOSS T ' ,IERR)
    elseif (tag(1:3) .eq. 'I23') then
       !   Sample input 
       !   '+I23 TN688 Far-Periphery Diffusion Rate ( < 0 = CDPERP ) '  -1.0
       !   IF (CDPERPFP.LT.0.0) CDPERPFP = CDPERP
       call divrd(CDPERPFP,.FALSE.,0.0,.FALSE.,0.0,'FP DIFF RATE  ' ,IERR)
       !  endif





    elseif     (tag(1:3) == 'I24') then 
       ! ----------------------------------------------------------------------- 

       !     TAG I24 : Ion initial positon in cell 

       !     Option to turn on/off more exact calculation of particle initial 
       !     positions. (default = 1 = on) 

       CALL ReadI(line,init_pos_opt,0,1,'Refined init position options') 

    else if (tag(1:3) == 'I25') then 
       ! ----------------------------------------------------------------------- 

       !     TAG I25 : FP Neutral ionization option 

       CALL ReadI(line,fp_neut_opt,0,1,'FP NEUTRAL IONIZATION OPTION') 

    else if (tag(1:3) == 'I26') then 
       ! ----------------------------------------------------------------------- 

       !     TAG I26 : FP plasma option 

       CALL ReadI(line,fp_plasma_opt,0,4,'FP PLASMA OPTION') 

    else if (tag(1:3) == 'I27') then 
       ! ----------------------------------------------------------------------- 

       !     TAG I27 : FP temperature 

       CALL ReadR(line,fp_te,0.0,HI,'FP Temperature in eV') 

    else if (tag(1:3) == 'I28') then 
       ! ----------------------------------------------------------------------- 

       !     TAG I28: FP density 

       CALL ReadR(line,fp_ne,0.0,HI,'FP Density in m-3') 

    else if (tag == 'I29') then 

       ! ----------------------------------------------------------------------- 



       !...  REPLACE! (?) - jdemod - not sure what is up here - why does it need 
       !                             replacement 

       !     TAG I29 and I30 - corner points for line/box injections 

       CALL Read2R(line,cxscA,cyscA,-HI,HI,'line/box injection point A') 

    else if (tag == 'I30') then 
       CALL Read2R(line,cxscB,cyscB,-HI,HI,'line/box injection point B') 

    else if (tag(1:3) == 'I31') then 
       ! ----------------------------------------------------------------------- 

       !     TAG I31: FP flow option 
       !     Far periphery transport flow option 
       !     0 - no flow in far periphery 
       !     1 - flow in far periphery is the same as associated ring 
       !     2 - flow in far periphery is specified as input 
       !     3 - flow pattern in far periphery is based on associated virtual ring 

       CALL ReadI(line,fp_flow_opt,0,3,'FP FLOW OPTION') 


    else if (tag(1:3) == 'I32') then 
       ! ----------------------------------------------------------------------- 

       !     TAG I32: FP flow velocity for option 2 

       CALL ReadR(line,fp_flow_velocity_input,-HI,HI,'FP Flow velocity for FP flow option 2') 

    else if (tag(1:3) == 'I33') then 
       ! ----------------------------------------------------------------------- 

       !     TAG I33: FP Number of radial bins in FP grid 

       CALL ReadI(line,fp_n_bins,1,maxnrs,'Number of radial bins for FP grid') 

    else if (tag(1:3) == 'I34') then 
       ! ----------------------------------------------------------------------- 

       !     TAG I34: FP Grid width option 

       !     Defines the option used to choose the width of simple crude 
       !     FP mesh 
       !     Option 0 = maximum distance from edge cell to wall from 
       !               fp_walldist values 
       !     Option 1 = width of grid is specified using fpxmaxo for MAIN and 
       !                fpxmaxi for PFZ 

       CALL ReadI(line,fp_grid_width_opt,0,1,'Number of radial bins for FP grid') 

    else if (tag(1:3) == 'I35') then 
       ! ----------------------------------------------------------------------- 

       !     TAG I35 : Main Chamber Ion reflection coefficient for fpopt 1 

       CALL ReadR(line,mc_recyc,0.0,1.0,'MC Reflection coefficient') 

    else if (tag(1:3) == 'I36') then 
       ! ----------------------------------------------------------------------- 

       !     TAG I36 : Private Fluz Zone Ion reflection coefficient for fpopt 1 

       CALL ReadR(line,pfz_recyc,0.0,1.0,'PFZ Reflection coefficient') 


    else if (tag(1:3) == 'I37') then 

       ! ----------------------------------------------------------------------- 

       !     TAG I37 : Fluid code charge state index to use for ion injection 
       !               options 12 and 13 

       CALL ReadI(line,e2diz_inj,1,maxe2dizs,'Fluid code impurity charge state index for injection') 

    else if (tag(1:3) == 'I38') then 
       ! Tag I38 for prompt deposition option 4, the average charge state 
       ! of the ion near the target. This value is used in the calculation 
       ! of the gyroradius instead of the calculated value. 
       call readr(line, prompt_dep_avg_z, 0.0, hi,'Average charge near target for prompt dep option 4') 

    ELSE 
       CALL ER('ReadTagSeriesI','Unrecognized tag',*99) 
    end if



    RETURN 
99  WRITE(SLOUT,'(5X,3A)') 'LINE = "',line,'"' 
    WRITE(SLOUT,'(5X,3A)') 'TAG  = "',tag ,'"' 
    WRITE(0    ,'(5X,3A)') 'LINE = "',line(1:LEN_TRIM(line)),'"' 
    WRITE(0    ,'(5X,3A)') 'TAG  = "',tag ,'"' 
    WRITE(0,*) '    DIVIMP HALTED' 
    STOP 
  END subroutine readtagseries_i

  ! ====================================================================== 




end module unstructured_input
