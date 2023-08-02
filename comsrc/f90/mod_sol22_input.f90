module mod_sol22_input


  implicit none
  !
  ! jdemod - I am creating this data module to start replacing common blocks.
  !
  ! Modules perform the same task but allow for both data sharing and initialization.
  ! In addition, executable code related to variable setup can also be included. 
  !
  ! Ultimately it would be good for all the unstructured input to be organized and added to files 
  ! like this one. 
  !

  integer,public :: debug_sol22 = 0        ! 284 - SOL22 debug flag - default 0 = off
  integer,public :: debug_sol22_ir = 1     ! 285 - SOL22 debug - ring for hi res plasma - default = 1
  integer,public :: debug_sol22_ikopt = 1  ! 286 - SOL22 debug - ikopt (ring end) for detailed plasma - default = 1

  

contains


  subroutine readsol(ierr)
    !use mod_params
    use mod_io
    use debug_options
    use mod_solparams
    use mod_solswitch
    use mod_solcommon

    implicit none

    !     This subroutine reads the input parameters for the case
    !     from the standard input or redirected from a file.

    !     include 'params'

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'

    integer ierr

    !      integer ierr

    !     Model parameters

    !real*8 spts(mxspts)

    call pr_trace('MOD_SOL22_INPUT','START OF READSOL')
    
    call rdi(forcet,.TRUE.,0,.true.,1,         'force te=ti'    ,ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER FORCET')

    CALL RDQ(initm0,.TRUE.,0.0d0,.FALSE.,0.0d0,'target mach num',IERR)
    CALL RDQ(deltam0,.TRUE.,0.0d0,.FALSE.,0.0d0,'delta mach num',IERR)

    !     ------------------------------------------------------------------

    !     Ionization Source

    CALL RDQ(m0res,.TRUE.,0.0d0,.FALSE.,0.0d0,'Resolution in m0',IERR)
    CALL RDI(lensind,.TRUE.,0,.TRUE.,1,     'ion source abs/rel',IERR)
    CALL RDQ(lensst,.TRUE.,0.0d0,.FALSE.,0.0d0,'ion src start',  IERR)
    CALL RDQ(lensfi,.TRUE.,0.0d0,.FALSE.,0.0d0,'ion src finish', IERR)
    CALL RDQ(lams,.TRUE. ,0.0d0,.FALSE.,0.0d0,'ion decay len   ',IERR)

    call pr_trace('MOD_SOL22_INPUT','AFTER IONIZATION SOURCE')
    !     Radiation source

    CALL RDQ(lenri,.TRUE.,0.0d0,.FALSE.,0.0d0,'rad source len  ',IERR)
    CALL RDQ(lamri,.TRUE. ,0.0d0,.FALSE.,0.0d0,'rad decay len  ',IERR)
    CALL RDQ(frri ,.TRUE. ,0.0d0,.FALSE.,0.0d0,'rad power mult ',IERR)
    CALL RDQ(alfimp,.TRUE.,0.0d0,.FALSE.,0.0d0,'nimp/ne ratio  ',IERR)
    CALL RDQ(talimp,.TRUE. ,0.0d0,.FALSE.,0.0d0,'base Temp     ',IERR)

    call pr_trace('MOD_SOL22_INPUT','AFTER RADIATION SOURCE')
    !     Miscellaneous

    CALL RDQ(ex1imp ,.FALSE. ,0.0d0,.FALSE.,0.0d0,'exponent 1  ',IERR)
    CALL RDQ(ex2imp ,.FALSE. ,0.0d0,.FALSE.,0.0d0,'exponent 2  ',IERR)
    CALL RDQ(gamcor,.false.,0.0d0,.FALSE.,0.0d0,'i power corr. ',IERR)

    CALL RDQ(gamecor,.false.,0.0d0,.FALSE.,0.0d0,'e power corr.',IERR)
    CALL RDQ(ceicf,.TRUE. ,0.0d0,.FALSE.,0.0d0, 'CX power frac ',IERR)

    CALL RDQ(recfrac,.TRUE.,0.0d0,.TRUE.,1.0d0,'Recycle frac ',IERR)
    call rdq(peicf,.true.  ,0.0d0,.false.,0.0d0,'Pei Correction',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER MISC')
    !     Power distribution

    call rdi(velsw,.true.  ,0,.true.,3,       'Vel Error Switch',ierr)
    call rdq(spowbeg,.true.,0.0d0,.true.,0.5d0,'Power Dist Beg',ierr)

    !     Gperp Distribution function

    call rdq(spowlen,.true.,0.0d0,.true.,0.5d0,'Power Dist Len',ierr)
    call rdq(gperpfrac,.true.,0.0d0,.true.,1.0d0,'Part Dist Frac',ierr)
    call rdq(gperpbegf,.true.,0.0d0,.true.,0.5d0,'Part Dist Beg',ierr)

    !     Extra Gperp source/sink - start and end positions.

    call rdq(gperpendf,.true.,0.0d0,.true.,0.5d0,'Part Dist Len',ierr)
    call rdq(gextra_mult,.true.,0.0d0,.false.,0.0d0,'Gextra flux mult',ierr)
    call rdq2(gextra_src_start,gextra_src_stop,.true.,0.0d0,.true.,1.0d0, 'Gextra SRC start/stop',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER GPERP')

    !     Field line length fraction for distributing the Private plasma
    !     electron and ion power loads.

    call rdq2(gextra_sink_start,gextra_sink_stop,.true.,0.0d0,.true.,1.0d0,'Gextra SINK start/stop',ierr)

    !     IK index for edge2d compatibility option 9

    call rdq(pp_pow_dist,.true.,0.0d0,.false.,0.0d0,'PP Pow Dist',ierr)
    call rdi(ike2d,.true.,1,.false.,0,'IK start Index for E2D-9',ierr)

    !     Cutoff temperature for PINQID term

    call rdi(fillopt,.true.,0,.true.,3,'Gap fill option - E2D-9',ierr)
    call rdq(tcutqe,.true.,0.0d0,.false.,0.0d0,'Cut T-PINQE',ierr)
    call rdq(tcutatiz,.true.,0.0d0,.false.,0.0d0,'Cut T-QIDATIZ',ierr)
    call rdq(tcutmliz,.true.,0.0d0,.false.,0.0d0,'Cut T-QIDMLIZ',ierr)
    call rdq(tcutrec,.true.,0.0d0,.false.,0.0d0,'Cut T- QIDREC',ierr)
    call rdq(tcutcx,.true.,0.0d0,.false.,0.0d0,'Cut T- QIDCX',ierr)
    call rdq(trefcx,.true.,0.0d0,.false.,0.0d0,'REF T- QIDCX 1',ierr)
    call rdq(tmin,.false.,0.0d0,.false.,0.0d0,'Min. Allowed T',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER PINQID')

    !     Momentum Source

    call rdq(dropfrac,.true.,0.0d0,.true.,1.0d0,'Allowed T-drop',ierr)
    call rdq(smom_mult,.false.,0.0d0,.false.,0.0d0,'Mom.Loss Multiplier',ierr)
    call rdq(ffric,.false.,0.0d0,.false.,0.0d0,'Mom.loss frac  ',ierr)
    CALL RDQ(lenmom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'mom source len',IERR)
    CALL RDQ(lammom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'mom decay len ',IERR)
    CALL RDQ(rcxmom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'cx/iz ratio   ',IERR)
    call rdq(tcxmom,.TRUE.,1.0001d0,.FALSE.,0.0d0,'T for CXmult',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER MOMENTUM')

    !     Source term multipliers

    call rdq(tcxcut,.true.,0.0d0,.false.,0.0d0,'Cut T -  CXmult',ierr)
    call rdq(qesrc_mult,.false.,0.0d0,.FALSE.,0.0d0,'PINQE mult',ierr)

    !     call rdq(qisrc_mult,.false.,0.0d0,.false.,0.0d0,'PINQI mult',ierr)

    !     ------------------------------------------------------------------

    call rdq(radsrc_mult,.false.,0.0d0,.FALSE.,0.0d0,'PINQE based PRAD mult',ierr)

    CALL RDI(ndiv,.TRUE., 1,.FALSE., 0,'NUMBER OF STEPS    ',IERR)

    call pr_trace('MOD_SOL22_INPUT','BEFORE SWITCHES')

    !     Read in switches 0.0 is off, 1.0 is on.

    call rdr(switch(swion),.true. ,0.0,.false.,0.0, 'alt ion',ierr)
    call rdr(switch(swioni),.true. ,0.0,.false.,0.0,'init ion',ierr)

    call rdr(switch(swionp),.false.,0.0,.false.,0.0,'pp ion',ierr)
    if (switch(swionp).eq.-1.0) then
       if (switch(swion).eq.1.0.or.switch(swion).eq.2.0.or.switch(swion).eq.8.0) then
          switch(swionp) = switch(swioni)
       else
          switch(swionp) = switch(swion)
       endif

    endif
    call rdr(switch(swcond),.true.,0.0,.false.,0.0, 'cond sw',ierr)
    call rdr(switch(swconv),.true.,0.0,.false.,0.0, 'conv sw',ierr)
    call rdr(switch(swprad),.true.,0.0,.false.,0.0, 'prad sw',ierr)
    call rdr(switch(swphelp),.true.,0.0,.false.,0.0,'phelp sw',ierr)
    call rdr(switch(swpei),.true.,0.0,.false.,0.0,  'pei sw ',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER BASIC SWITCHES')

    !     PINQID switches

    call rdr(switch(swpcx),.true.,0.0,.false.,0.0,  'pcx sw ',ierr)
    call rdr(switch(swqidatiz),.true.,0.0,.false.,0.0,'atiz sw',ierr)
    call rdr(switch(swqidmliz),.true.,0.0,.false.,0.0,'mliz sw',ierr)
    call rdr(switch(swqidrec),.true.,0.0,.false.,0.0,'rec sw',ierr)

    call rdr(switch(swqidcx),.true.,0.0,.false.,0.0,'cx sw',ierr)
    call rdr(switch(swppelec),.true.,0.0,.false.,0.0,'pp elec sw',ierr)


    !     jdemod -
    !     switch(swppress) is read in using tag 283 of unstructured input
    !     default value is 0 or OFF

    !      call rdr(switch(swppress),.true.,0.0,.false.,0.0,'pp power sw',ierr)


    call rdr(switch(swppion),.true.,0.0,.false.,0.0,'pp ion sw',ierr)
    call rdr(switch(swvisc1),.true.,0.0,.true.,0.0,'visc1 sw',ierr)
    call rdr(switch(swnmom),.true.,0.0,.false.,0.0, 'N mom sw',ierr)

    call pr_trace('MOD_SOL22_INPUT','HALFWAY THROUGH SWITCHES')


    !     Read in Edge2d compatibility switch and the subsequent values of
    !     ne, Te and Ti at the centre point of the first cell

    call rdr(switch(swmach),.true.,0.0,.false.,0.0, 'mach sw',ierr)
    call rdr(switch(swe2d),.false.,0.0,.false.,0.0, 'e2d sw',ierr)
    call rdr(switch(swpow),.true.,0.0,.false.,0.0, 'power sw',ierr)

    call rdr(switch(swpowp),.false.,0.0,.false.,0.0, 'pp pow',ierr)

    if (switch(swpowp).eq.-1.0) switch(swpowp) = switch(swpow)
    call rdr(switch(swgperp),.true.,0.0,.false.,0.0,'GamPerp',ierr)

    call rdr(switch(swgperpp),.false.,0.0,.false.,0.0,'GamPerpP',ierr)

    !     Extra Gperp source/sink term

    if (switch(swgperpp).eq.-1.0) switch(swgperpp) = switch(swgperp)

    call rdr(switch(swextra),.true.,0.0,.false.,0.0,'GP Src/Sink',ierr)
    call rdr(switch(swmajr),.true.,0.0,.false.,0.0,'MajorRad',ierr)
    call rdr(switch(swcore),.true.,0.0,.false.,0.0,'Core Src',ierr)
    call rdr(switch(swrecom),.true.,0.0,.false.,0.0,'Recomb.',ierr)
    call rdr(switch(swsmooth),.true.,0.0,.false.,0.0,'Smooth',ierr)
    call rdr(switch(swdetach),.true.,0.0,.false.,0.0,'Detach',ierr)
    call rdr(switch(swerror),.true.,0.0,.false.,0.0,'ERROR',ierr)

    call pr_trace('MOD_SOL22_INPUT','DONE SWITCHES')


    !CALL RDRARN(deflist,ndef,mxspts,0.0,real(maxnrs),.FALSE.,0.0,MACHHI,2,'DEFAULT SOLVER DATA',IERR)
    !
    ! jdemod - The max value should the maximum number of rings but to avoid dependency on params
    !          module I am removing this constraint. 
    CALL RDRARN(deflist,ndef,mxspts,0.0,SOL22_HI,.FALSE.,0.0,SOL22_MACHHI,2,'DEFAULT SOLVER DATA',IERR)

    !     Default plots to off - this can be changed in the calcsol_interface
    !     routine in solascv1.f

    graph = 0
    graphaux = 0
    graphvel = 0

    !      CALL RDI(graph,.TRUE., 0,.TRUE., 1,'GRAPH OPTION       ',IERR)
    !      CALL RDI(graphaux,.TRUE.,0,.TRUE.,1,'AUX GRAPH OPTION  ',IERR)
    !      CALL RDI(graphvel,.TRUE.,0,.TRUE.,1,'VEL GRAPH OPTION  ',IERR)
    !      call rdr(graphran,.true.,0.0,.false.,0.0,'CXMAX VALUE  ',ierr)


    graphran = 0.0

    call pr_trace('MOD_SOL22_INPUT','END OF READSOL')


    return



  end subroutine readsol

  subroutine sol22_initialize_unstructured_input
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none

    
!
!      TAG 282: SOL22
!
!     Initialization of Array input for tag 282 specifying ffric 
!     values on a ring by ring basis for both targets. 
!
      n_extffric = 0
!      extffric = 0.0
!     
!     TAG 283: SOL22 - private plasma pressure loss option
!
!     Set the default for this value to OFF = 0
!
      switch(swppress) = 0.0
!     
!     TAG 284: SOL22 - debug SOL22 
!
!     Set the default for this value to OFF = 0
!
      debug_sol22 = 0
!     
!     TAG 285: SOL22 - debug SOL22 
!
!     Set the default for this value to OFF = 0
!
      debug_sol22_ir = 1
!     
!     TAG 286: SOL22 - debug SOL22 
!
!     Set the default for this value to OFF = 0
!
      debug_sol22_ikopt = 1
!
!     TAG 287: SOL22 - base ionization source length for algorithmic ionization options
!
      alg_ion_src_len = 2.0
!
!     TAG 288: SOL22 - ring by ring specification of radiation loss parameters
!
!     Initialization of Array input for tag 288 specifying radiation parameters
!     values on a ring by ring basis for both targets. 
!
      n_extradsrc = 0
      extradsrc = 0.0 
      !
      ! TAG 289: swqperpe - electron perpendicular power flux option
      !          compensates for power terms giving conservation of energy on flux tube
      !          default = OFF
      !
      switch(swqperpe) = 0.0

      !
      ! TAG 290: swqperpi - electron perpendicular power flux option
      !          compensates for power terms giving conservation of energy on flux tube
      !          default = OFF
      !
      switch(swqperpi) = 0.0

      !
      ! TAG 291: DBLSRC_OPT - select option for double source 0=exp+rect  1=rect+rect
      !
      dblsrc_opt = 0

      !
      ! TAG 292: DBLSRC_FRAC - fraction in first source - rest is 1-fraction
      !
      dblsrc_frac = 0.5

      !
      ! TAG 293: DBLSRC 1 - Opt 0 = lams len   Opt 1 = len1 len2
      !
      dblsrc1_p1 = 0.05
      dblsrc1_p2 = 0.1
      !
      ! TAG 294: DBLSRC 2   len1 len2
      !
      dblsrc2_p1 = 0.0
      dblsrc2_p2 = 0.5

      !
      ! TAG 295: switch(swepow) - external electron power term
      !
      switch(swepow) = 0.0

      !
      ! TAG 296: switch(swipow) - external ion power term
      !
      switch(swipow) = 0.0

      !
      ! TAG 297: switch(swepow)=2.0 - external electron power term file name
      !
      ext_epow_fn = 'ext_epow_data.txt'

      !
      ! TAG 298: switch(swipow)=2.0 - external ion power term file name
      !
      ext_ipow_fn = 'ext_ipow_data.txt'

      !
      ! TAG 299: Sol22 print option - used to select more detailed sol22 outputs in the .lim file (fort.6)
      !    sol22_print =0 off
      !
      sol22_print = 0
      
  end subroutine sol22_initialize_unstructured_input


  subroutine sol22_unstructured_input(tag,line,ierr)
    use mod_io
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    character*(*) :: tag,line
    integer :: ierr

    !
    !
    !
    ! -----------------------------------------------------------------------
    !
    !     Options added by David for various parts of the code. 
    !
    ! -----------------------------------------------------------------------
    !
    !
    ! -----------------------------------------------------------------------
    !
    !     Trying to keep options in alphabetical order. 
    !
    if (tag(1:3).eq.'282') then  
    !     TAG 282 - SOL option 22 - reads in FFRIC values for 
    !               the momentum loss option so that they may vary
    !               from ring to ring and target to target.
    !             - the format is IR FFRIC1 FFRIC2
    !               jdemod - format is IR FFRIC1 LEN1 FFRIC2 LEN2 for option 1
    !                      - format is IR FFRIC1 LAM1 FFRIC2 LAM2 for option 2 - LEN is set to SMAX/2
    !
    !               FFRIC1 etc applies to the first half ring - this would
    !               be the OUTER half ring for X-point up grids and the 
    !               INNER half for X-point down grids. This may also
    !               be designated target 2. 
    !               FFRIC2 etc applies to the second half of the ring  
    !     
    !     Note: the tag line precedes a standard DIVIMP array input of
    !           three lines.  
    !
       CALL RDQARN(extffric,n_extffric,MXSPTS,-sol22_MACHHI,sol22_MACHHI,.FALSE.,&
               -sol22_machhi,sol22_MACHHI,4,'SET OF MOM-LOSS COEF BY RING',IERR)
          
    elseif (tag(1:3).eq.'283') then  
       !
       !     jdemod
       !     TAG 283 - SOL option 22 - reads the value for switch(swppress) 
       !             - this is the pfz target pressure redistribution option. 
       !     

       CALL ReadR(line,switch(swppress),0.0,2.0,'SOL22 PFZ PRESSURE LOSS OPT')
       !
       !        write(0,*) '*283 - Read switch',swppress,switch(swppress)
       !
    elseif (tag(1:3).eq.'284') then  
       !     jdemod
       !     TAG 284 - SOL option 22 - SOL22 debugging switch
       !     

       CALL ReadI(line,debug_sol22,0,1,'SOL22: DEBUG SWITCH')
       !
    elseif (tag(1:3).eq.'285') then  
       !     jdemod
       !     TAG 285 - SOL option 22 - SOL22 debug IR for detailed profile
       !     

       ! jdemod
       ! maxnrs might not be defined depending on the usage so this is just an index.
       ! Make imax very large for now
       !CALL ReadI(line,debug_sol22_ir,1,maxnrs,'SOL22 DEBUG RING')
       CALL ReadI(line,debug_sol22_ir,1,100000,'SOL22 DEBUG RING')
       !
    elseif (tag(1:3).eq.'286') then  
       !
       !     jdemod
       !     TAG 286 - SOL option 22 - SOL22 debug IKOPT for detailed profile
       !     

       CALL ReadI(line,debug_sol22_ikopt,1,2,'SOL22 DEBUG IKOPT-RING END')
    elseif (tag(1:3).eq.'287') then  
       !
       !     TAG 287: SOL22 - base ionization source length for algorithmic ionization options
       !

       CALL ReadR(line,alg_ion_src_len,0.0,sol22_MACHHI,'DEFAULT IONIZATION SOURCE LENGTH FOR ALGORITHMIC OPTIONS')
    elseif (tag(1:3).eq.'288') then  
       !
       !     TAG 288 - SOL option 22 - reads in radiation values for 
       !               the radiation option so that they may vary
       !               from ring to ring and target to target.
       !             - the format is IR LENR1 LAMR1 FFR1 LENR2 LAMR2 FFR2
       !
       !               "1" applies to the first half ring - this would
       !               be the OUTER half ring for X-point up grids and the 
       !               INNER half for X-point down grids. This may also
       !               be designated target 2. 
       !               "2" applies to the second half of the ring  
       !     
       !     Note: the tag line precedes a standard DIVIMP array input of
       !           three lines.  
       !
       CALL RDQARN(extradsrc,n_extradsrc,MXSPTS,-sol22_MACHHI,sol22_MACHHI,.FALSE.,&
            -sol22_machhi,sol22_MACHHI,6,'SET OF RADIATION COEF BY RING',IERR)
       !
    elseif (tag(1:3).eq.'289') then  
       !
       !     jdemod
       !     TAG 289 - SOL option 22 - reads the value for switch(swqperpe) 
       !             - switch to balance electron power on flux tube
       !             - default is OFF
       !     

       CALL ReadR(line,switch(swqperpe),0.0,1.0,'SOL22: ELECTRON PERPENDICULAR POWER FLUX OPTION')
       !
    elseif (tag(1:3).eq.'290') then  
       !
       !     jdemod
       !     TAG 290 - SOL option 22 - reads the value for switch(swqperpi) 
       !             - switch to balance ion power on flux tube
       !             - default is OFF
       !     

       CALL ReadR(line,switch(swqperpi),0.0,1.0,'SOL22: ION PERPENDICULAR POWER FLUX OPTION')
       !
    elseif (tag(1:3).eq.'291') then  
       !
       !     jdemod  tag 291-294
       !     Double shape ionization source
       !     Option 0 : exp + rect (default)
       !     Option 1 : rect + rect
       !
       !     Other paramters
       !     dblsrc_frac : fraction of particles in first source
       !     dblsrc1_p1, dlbsrc1_p2 : paramters for first source - either lams len or len1 len2
       !     dblsrc2_p1, dlbsrc2_p2 : paramters for second source - len1 len2
       !     
       CALL ReadI(line,dblsrc_opt,0,1,'OPTION for dual profile ionization source')

       !
    elseif (tag(1:3).eq.'292') then  
       !
       !     jdemod
       !     TAG 292
       !     dblsrc_frac : fraction of particles in first source
       !     
       CALL ReadR(line,dblsrc_frac,0.0,1.0,'SOL22: DBLSRC_FRAC : fraction of particles in first source')
       !
    elseif (tag(1:3).eq.'293') then  
       !
       !     jdemod
       !     TAG 293
       !     dblsrc1_p1, dlbsrc1_p2 : paramters for first source - either lams len or len1 len2
       !     
       CALL Read2R(line,dblsrc1_p1,dblsrc1_p2,0.0,1.0,'SOL22: DBLSRC1_P1,2 : parameters for first source')
       !
    elseif (tag(1:3).eq.'294') then  
       !
       !     jdemod
       !     TAG 294
       !     dblsrc2_p1, dlbsrc2_p2 : paramters for second source len1 len2
       !     
       CALL Read2R(line,dblsrc2_p1,dblsrc2_p2,0.0,1.0,'SOL22: DBLSRC2_P1,2 : parameters for second source')
       !
    elseif (tag(1:3).eq.'295') then  
       !
       !     jdemod
       !     TAG 295 - SOL option 22 - reads the value for switch(swepow) for external electron power term 
       !             - default is OFF = 0.0
       !             - option 1 = read data from div aux input file - Looks for tag EXTEPOW:
       !                          Add dataset to the divimp aux input file using the following format
       ! 
       !                          EXTEPOW:
       !                          format(6e18.10)   ((data_array(ik,ir,1),ik=1,nks(ir)),ir=1,nrs)
       ! 
       !             - option 2 = read raw data from R,Z,EPOW file and interpolate onto DIVIMP mesh
       !                          (external data source file must be connected to ext_epow_data.txt)
       !                          External file format is:
       !                          # or $ = Comments ignored on file read
       !                          NROWS: <nr>           
       !                          NCOLS: <nz> 
       !                          DATA:
       !                           R Z RAD       These may be space or comma separated
       !               Any grid coordinates outside the supplied data will be assigned a zero           
       !               Data must be entered from low R to high R and low Z to high Z
       !               i.e. Start from lower left corner of data set
       !               Code could be added to support additional data formats or organization 
       !     

       CALL ReadR(line,switch(swepow),0.0,2.0,'SOL22: EXTERNAL ELECTRON POWER TERM OPTION')
       !
    elseif (tag(1:3).eq.'296') then  
       !
       !     jdemod
       !     TAG 295 - SOL option 22 - reads the value for switch(swipow) for external ion power term 
       !             - default is OFF = 0.0
       !             - option 1 = read data from div aux input file - Looks for tag EXTIPOW:
       ! 
       !                          EXTIPOW:
       !                          format(6e18.10)   ((data_array(ik,ir,1),ik=1,nks(ir)),ir=1,nrs)
       ! 
       !             - option 2 = read raw data from R,Z,IPOW file and interpolate onto DIVIMP mesh
       !                          (external data source file must be connected to the ext_ipow_data.txt)
       !                          External file format is:
       !                          # or $ = Comments ignored on file read
       !                          NROWS: <nr>           
       !                          NCOLS: <nz> 
       !                          DATA:
       !                           R Z RAD       These may be space or comma separated
       !               Any grid coordinates outside the supplied data will be assigned a zero           
       !               Data must be entered from low R to high R and low Z to high Z
       !               i.e. Start from lower left corner of data set
       !               Code could be added to support additional data formats or organization 
       !     

       CALL ReadR(line,switch(swipow),0.0,2.0,'SOL22: EXTERNAL ION POWER TERM OPTION')
       !
    elseif (tag(1:3).eq.'297') then  
       !
       !     jdemod
       !     TAG 297 - SOL option 22 - file name with external electron power source for swepow option = 2.0
       !                               default = ext_epow_data.txt 
       !     

       CALL ReadC(line,ext_epow_fn,'SOL22: EXTERNAL ELECTRON POWER TERM FILE NAME')

       !
    elseif (tag(1:3).eq.'298') then  
       !
       !     jdemod
       !     TAG 298 - SOL option 22 - file name with external ion power source for swipow option = 2.0
       !                               default = ext_ipow_data.txt 
       !     

       CALL ReadC(line,ext_ipow_fn,'SOL22: EXTERNAL ION POWER TERM FILE NAME')

       !

    elseif (tag(1:3).eq.'299') then  

       ! jdemod
       ! TAG 299: Sol22 print option - used to select more detailed sol22 outputs in the .lim file (fort.6)
       !    sol22_print =0 off
       !
       CALL ReadI(line,sol22_print,0,2,'Selects level of detail in SOL22 output to the .lim(fort.6) file ')

    endif

  end subroutine sol22_unstructured_input

end module mod_sol22_input
