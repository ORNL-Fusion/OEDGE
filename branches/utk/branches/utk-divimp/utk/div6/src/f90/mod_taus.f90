module taus

  !
  ! Mod_taus: This module replaces the characteristic time calculation code 
  !           from the tau.f source code file. This code is 
  !           self-contained after the init call is made setting the required
  !           parameters. 
  !
  ! jde - June 20, 2005
  !

  !
  ! retain local variable values across multiple uses of the module 
  !
  save

  real,private ::      c215a,c350a,c350b
  REAL,private ::      ROOTMI,ROOTTT
  real,private,PARAMETER :: LAMBDA=15.0
  REAL,private ::      FTAU,FTAUP,FTAUS,FTAUT,RIZSQR,STAU,TAU

  real,private ::      iz_bg_plasma,zeff

  integer,private :: coll_opt     ! cioptb
  integer,private :: stop_opt     ! cioptc
  integer,private :: heat_opt     ! cioptd
  integer,private :: spec_ring    ! irspec

  real,private :: mb,mi           ! crmb, crmi (or crmhc)
  real,private :: ave_temp        ! ctemav  

  real,private :: timestep        ! qtim

  real,private :: last_kfps,last_kkkfps,last_kfss,last_kfts  ! record values set on last call to eval_taus

  contains

  subroutine init_taus(crmb,crmi,rizb,cioptb,cioptc,cioptd,czenh,cizeff,ctemav,irspec,qtim)
    implicit none
    real czenh,crmb,crmi,rizb,qtim,ctemav
    integer cioptb,cioptc,cioptd,irspec,cizeff

    !
    ! Calculate base expressions used in the evaluation of the transport coefficients
    !

    FTAU  = CZENH * SQRT(CRMB) * rizb * rizb * LAMBDA * QTIM
    FTAUP = FTAU * 6.8E-14
    FTAUS = FTAU * 6.8E-14 * (1.0 + CRMB/CRMI)
    FTAUT = FTAU * 1.4E-13

    C215A = FTAU * 1.4E-13 * SQRT(CRMI)
    C350A = 9.0E-14 * (1.0+CRMI/CRMB) * rizb * rizb  * CZENH * LAMBDA * QTIM / SQRT(CRMI)
    C350B = 9.0E-14 * SQRT(CRMI) * rizb * rizb * CZENH * LAMBDA * QTIM /  CRMB
    

    iz_bg_plasma = rizb
    zeff = real(cizeff)
    ave_temp=ctemav
    roottt = sqrt(ctemav)

    !
    ! Copy options
    ! 
    coll_opt = cioptb
    stop_opt= cioptc
    heat_opt = cioptd
    spec_ring = irspec

    !
    ! Copy mass and timestep 
    !
    mb = crmb
    call taus_set_mass(crmi)

    timestep = qtim


  end subroutine init_taus


  subroutine taus_set_mass(massi)
    implicit none
    real massi

    mi = massi
    ROOTMI = SQRT (massi)

  end subroutine taus_set_mass



  subroutine eval_taus(ik,ir,iz,ne,ti,kfps,kkkfps,kfss,kfts)

    use mod_params

    implicit none
    integer ik,ir,iz
    real ne,ti,kfps,kkkfps,kfss,kfts

    ! local variables
    real rizsqr 
    real stau
    


    RIZSQR = REAL (IZ) * REAL (IZ)

    if (ti.ne.0.0) then 
       STAU = ne / (MI * ti**1.5) *RIZSQR
    else
       write(0,*) 'WARNING: EVAL_TAUS: Ti=0.0'
       stau = 1.0e10
    endif

    !
    !-----------------------------------------------------------------------
    !           TAU PARALLEL           NOTES 3,50,103
    !-----------------------------------------------------------------------
    !
    !---------- STANDARD CASE:-
    !---------- CFPS = 2.DELTAT.TI/TAUPARALLEL.
    !---------- FOR NON-ZERO OPTIONS, SPECIAL CASE APPLIES FOR K > KSPEC
    !---------- ONLY.  EACH TIME AN ION ENTERS THIS REGION ITS
    !---------- TEMPERATURE IS COMPARED AGAINST A PREVIOUS VALUE TO SEE
    !---------- WHETHER ANY VALUES NEED RECALCULATING.
    !

    KFPS = STAU * ti * FTAUP * 2.0

    if (ave_temp.eq.0.0.or.roottt.eq.0.0.or.rootmi.eq.0.0.or.mb.eq.0.0.or.ti.eq.0.0.or.mi.eq.0.0) then 
       write(0,'(a,6(1xg18.8))') 'WARNING:EVAL_TAUS: UNEXPECTED 0.0 VALUE:', roottt,rootmi,ave_temp,mb,ti,mi
    endif
    
    IF     (COLL_OPT.EQ.1.AND.IR.GE.SPEC_RING) THEN
       KFPS = 0.0
    ELSEIF (COLL_OPT.EQ.2.AND.IR.GE.SPEC_RING) THEN
       KFPS = 2.0 * ne * 6.8E-14 *  IZ_BG_PLASMA * ZEFF * RIZSQR * LAMBDA/(ROOTMI * ROOTTT) * TIMESTEP
    ELSEIF (COLL_OPT.EQ.4.AND.IR.GE.SPEC_RING.AND.AVE_TEMP.GT. TI *MI/MB) THEN
       KFPS = C350B * RIZSQR * ti *   ne / (AVE_TEMP * ROOTTT)
    ENDIF


    !
    !---------- SET ADDITIONAL ARRAY FOR USE IN INNERMOST LOOP OF LIM2
    !---------- EQUIVALENT TO SOME CONSTANT TIMES CFPS VALUES
    !---------- THIS SAVES 20% OF CPU TIME BY ELIMINATING A SQUARE ROOT
    !---------- KKKFPS = SQRT (4.88E8/(KFPS.MI)) . DELTAT. SIN(THETAB)
    !---------- FOR NOTE 284, SET KKKFPS = AS ABOVE . SQRT(2TI/KFPS)
    !

    IF (KFPS.EQ.0.0) THEN
       KKKFPS = 0.0
    ELSEIF (  (COLL_OPT.EQ.3.or.coll_opt.eq.7).AND.IR.GE.SPEC_RING) THEN
       KKKFPS = SQRT (9.76E8 * AVE_TEMP / MI) *  TIMESTEP / KFPS
    elseif (coll_opt.eq.5) then
       !
       !              Base diffusive velocity step size.
       !
       KKKFPS = 1.56e4 * SQRT(AVE_TEMP/MI) * TIMESTEP
       !
    elseif (coll_opt.eq.6) then
       !
       !              The velocity step for collision option 6 becomes
       !              (8kT/PIMi)**1/2 * (deltat/taupara)**1/2
       !              However taupara is proportional to Ti so the Ti in 8kTi
       !              and the Ti in taupara cancel leaving deltaV a constant
       !              for collision option 6. This saves CPU time.
       !              The effect is similar for collision option 11.
       !
       KKKFPS = 1.56e4 *  SQRT(1.0/MI *kfps/2.0) * TIMESTEP
       !
       !
    elseif (coll_opt.eq.8) then
       !
       !              Equivalent to collision option 0.
       !
       !              Correction factors have been added shifting from
       !              8kT/PIm to kT/m and also correcting for a factor
       !              of sqrt(2.0) ... this results in a change of
       !              PI/4.0 inside the square root.
       !
       KKKFPS = TIMESTEP * SQRT ((PI/4.0 * 4.88E8) /(KFPS*MI))
       !
    elseif (coll_opt.eq.9.and.ir.ge.spec_ring) then
       !
       !              Equivalent to collision option 3
       !
       KKKFPS=SQRT(((PI/4.0)* 9.76E8)*AVE_TEMP/MI)*TIMESTEP / KFPS
       !
    elseif (coll_opt.eq.10) then
       !
       !              Equivalent to collision option 5
       !
       !
       !              Base diffusive velocity step size.
       !
       KKKFPS = 1.56e4 * SQRT(PI/4.0 * AVE_TEMP/MI)   * TIMESTEP
    elseif (coll_opt.eq.11.or.coll_opt.eq.12) then
       !
       !              Equivalent to collision option 6
       !
       !
       !              The velocity step for collision option 11 becomes
       !              (kT/Mi)**1/2 * ( 2 deltat/taupara)**1/2
       !              However taupara is proportional to Ti so the Ti in kTi
       !              and the Ti in taupara cancel leaving deltaV a constant
       !              for collision option 11. This saves CPU time.
       !              The effect is similar for collision option 6. This
       !              differs from collision option 6 in the numerical
       !              constant sqrt(PI/4.0).
       !
       KKKFPS = 1.56e4 *SQRT(PI/4.0 * 1.0/MI *kfps/2.0) * TIMESTEP
       !
       !
    elseif (coll_opt.eq.13) then
       !
       !              The velocity step for collision option 13 becomes
       !              (kT/Mi)**1/2 * ( 2 deltat/taupara)**1/2
       !              However taupara is proportional to Ti so the Ti in kTi
       !              and the Ti in taupara cancel leaving deltaV a constant
       !              for collision option 11. This saves CPU time.
       !              The effect is similar for collision option 6. This
       !              differs from collision option 6 in the numerical
       !              constant sqrt(PI/4.0).
       !
       !              In option 13 ... Taupara is divided by a correction
       !              factor of (1+MB/MI).
       !
       !              This means kfps * (1+Mb/Mi)
       !
       KKKFPS = 1.56e4 * SQRT(PI/4.0 * 1.0/MI *(kfps*(1.0+MB/MI))  /2.0) * TIMESTEP
       !
       !
    ELSE
       KKKFPS = TIMESTEP * SQRT (4.88E8 /(KFPS*MI))
    ENDIF

    !
    !-----------------------------------------------------------------------
    !           TAU STOPPING           NOTES 3,50,103
    !-----------------------------------------------------------------------
    !
    !---------- STANDARD CASE:-
    !---------- KFSS = 1 - EXP (-DELTAT/TAUSTOPPING)
    !---------- NON ZERO OPTIONS APPLY FOR K > KSPEC ONLY ...
    !---------- FOR OPTION 2, TAUSTOPPING = TAUPARALLEL
    !----------                           = 2.DELTAT.TI/CFPS
    !

    TAU = STAU * FTAUS
    IF (TAU.GT.1.E-3) THEN
       KFSS = 1.0 - EXP(-TAU)
    ELSE
       KFSS = TAU
    ENDIF
    
    IF     (STOP_OPT.EQ.1.AND.IR.GE.SPEC_RING) THEN
       KFSS = 0.0
    ELSEIF (STOP_OPT.EQ.2.AND.IR.GE.SPEC_RING) THEN
       TAU = KFPS / (2.0*AVE_TEMP)
       IF (TAU.GT.1.E-3) THEN
          KFSS = 1.0 - EXP(-TAU)
       ELSE
          KFSS = TAU
       ENDIF
    ELSEIF (STOP_OPT.EQ.3.AND.IR.GE.SPEC_RING.AND.AVE_TEMP.GT. TI *MI/MB) THEN
       TAU = C350A * RIZSQR *  NE  / (AVE_TEMP * ROOTTT)
       IF (TAU.GT.1.E-3) THEN
          KFSS = 1.0 - EXP(-TAU)
       ELSE
          KFSS = TAU
       ENDIF
    ENDIF
    !
    !-----------------------------------------------------------------------
    !           TAU HEATING            NOTES 3,50,103,215
    !-----------------------------------------------------------------------
    !
    !---------- STANDARD CASE:-
    !---------- KFTS = 1 - EXP (-DELTAT/TAUHEATING)
    !---------- NON ZERO OPTIONS APPLY FOR K > KSPEC ONLY ...
    !
    TAU = STAU * FTAUT
    IF (TAU.GT.1.E-3) THEN
       KFTS = 1.0 - EXP(-TAU)
    ELSE
       KFTS = TAU
    ENDIF

    IF     (HEAT_OPT.EQ.1.AND.IR.GE.SPEC_RING) THEN
       KFTS = 0.0
    ELSEIF (HEAT_OPT.EQ.2.AND.IR.GE.SPEC_RING) THEN
       KFTS = 1.0
    ELSEIF (HEAT_OPT.EQ.3) THEN
       TAU = C215A * RIZSQR *  NE  /   ((MI *  TI  + MB * AVE_TEMP) ** 1.5)
       IF (TAU.GT.1.E-3) THEN
          KFTS = 1.0 - EXP(-TAU)
       ELSE
          KFTS = TAU
       ENDIF
    ENDIF

    last_kfps   = kfps
    last_kkkfps = kkkfps
    last_kfss   = kfss
    last_kfts   = kfts

  end subroutine eval_taus


  subroutine adjust_taus(ik,ir,iz,temi,lfps,lllfps,lfss,lfts,ne,tib)
  
      use hc_get
      use hc_put

      implicit none
      real :: temi  ! ion temperature
      real :: lfps,lllfps,lfss,lfts,ne,tib
      integer ik,ir,iz
      ! 
      ! Local variables
      !
      !real,external :: gltolds
      real :: temold,tau
      real :: ratio1,ratio2


!
!-----------------------------------------------------------------------
!       PARALLEL DIFFUSION AND ION TEMPERATURE
!-----------------------------------------------------------------------
!
!------ QUICKLY ADJUST RELEVANT CHARACTERISTIC TIMES
!------ FOR NON-STANDARD PLASMA OPTIONS.  BECAUSE OF THE SQUARE ROOTS
!------ THESE CALCULATIONS ARE DONE SPARINGLY TO PREVENT DRAMATIC
!------ SLOWING OF PROGRAM EXECUTION - FOR EXAMPLE, WHENEVER THE
!------ TEMPERATURE CHANGES BY 20% OR MORE.
!
        IF( ((COLL_OPT.eq.2.OR.coll_opt.eq.3.or.coll_opt.eq.4.or.coll_opt.eq.7.or.coll_opt.eq.9.or.&
            & STOP_OPT.eq.2.or.stop_opt.eq.3).AND.IR.GE.SPEC_RING).OR.HEAT_OPT.GE.3.or.coll_opt.eq.5.or.coll_opt.eq.10) THEN

            TEMOLD = GLTOLDS(IK,IR,IZ)
            IF (TEMI.GT.1.2*TEMOLD .OR. TEMI.LT.0.8*TEMOLD) THEN

              call PLTOLDS(IK,IR,IZ, TEMI)
!
              IF (COLL_OPT.EQ.2) THEN
                RATIO1 = SQRT (TEMOLD / TEMI)
                RATIO2 = 1.0 / SQRT (RATIO1)

                LFPS = RATIO1 * LFPS
                LLLFPS = RATIO2 * LLLFPS

              ELSEIF (COLL_OPT.EQ.3.or.coll_opt.eq.7.or.coll_opt.eq.9.or.coll_opt.eq.5.or.coll_opt.eq.10) THEN
                RATIO2 = SQRT (TEMI / TEMOLD)
                LLLFPS = RATIO2 * LLLFPS
              ELSEIF (COLL_OPT.EQ.4) THEN
                IF (TEMI.LE.TIB*MI/MB) THEN
                  LFPS = last_KFPS
                  LLLFPS = last_KKKFPS
                ELSE
                  LFPS = C350B*REAL(IZ*IZ)*TIB *  ne / (TEMI ** 1.5)
                  LLLFPS = TIMESTEP * SQRT(4.88E8/(LFPS*MI))
                ENDIF
              ENDIF
!
              IF (STOP_OPT.EQ.2) THEN
                TAU = LFPS / (2.0*TEMI)
                IF (TAU.GT.1.E-3) THEN
                  LFSS = 1.0 - EXP(-TAU)
                ELSE
                  LFSS = TAU
                ENDIF
              ELSEIF (STOP_OPT.EQ.3) THEN
                IF (TEMI.LE.TIB*MI/MB) THEN
                  LFSS = last_KFSS
                ELSE
                  TAU = C350A * REAL(IZ*IZ) * NE / (TEMI ** 1.5)
                  IF (TAU.GT.1.E-3) THEN
                    LFSS = 1.0 - EXP(-TAU)
                  ELSE
                    LFSS = TAU
                  ENDIF
                ENDIF
              ENDIF
!
              IF (HEAT_OPT.EQ.3) THEN
                TAU = C215A * REAL(IZ*IZ) * NE / ((MI * TIB + MB * TEMI) ** 1.5)
                IF (TAU.GT.1.E-3) THEN
                  LFTS = 1.0 - EXP (-TAU)
                ELSE
                  LFTS = TAU
                ENDIF
              ENDIF
            ENDIF
        ENDIF


  end subroutine adjust_taus



end module taus
