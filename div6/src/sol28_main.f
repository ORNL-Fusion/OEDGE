c     -*-Fortran-*-
c
c ====================================================================
c
      SUBROUTINE AssignCorePlasma(itube,nnode,mnode,node)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      INTEGER, INTENT(IN) :: itube,nnode,mnode
      TYPE(type_node)     :: node(nnode)

      
      INTEGER ion,ic1,ic2

      IF (nnode.NE.3.OR.mnode.NE.2) 
     .  CALL ER('AssignCorePlasma','Non-standard core node setup',*99)

      ion = 1

      ic1 = tube(itube)%cell_index(LO)
      ic2 = tube(itube)%cell_index(HI)

      fluid(ic1:ic2,ion)%ne = node(mnode)%ne
      fluid(ic1:ic2,ion)%te = node(mnode)%te
      fluid(ic1:ic2,ion)%ni = node(mnode)%ne
      fluid(ic1:ic2,ion)%ti = node(mnode)%ti(ion)

      RETURN
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE ProcessKineticIons
      USE mod_sol28
      USE mod_sol28_solver
      IMPLICIT none


      SELECTCASE (0)
        CASE (0)
c...      None at the moment:
          ni(0:icmax+1,0) = 0.0
        CASE DEFAULT
      ENDSELECT


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE PrescribeTemperatureProfiles
      USE mod_sol28
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER ion,ic1,ic2,target
      LOGICAL cont

      cont = .TRUE.
 
      DO WHILE(cont)
        cont = .FALSE.

c...    Te:
        CALL InterpolateTeProfile(1)

c...    Ti:
        ion = 1
        DO target = LO, HI
          ic1 = icbnd1(target)
          ic2 = icbnd2(target)
          IF (target.EQ.LO.AND.ic1.EQ.1    ) ic1 = 0
          IF (target.EQ.HI.AND.ic2.EQ.icmax) ic2 = icmax + 1
          SELECTCASE (opt%ti(target))
            CASE (0)
              ti(ic1:ic2,ion) = te(ic1:ic2) * DBLE(opt%ti_ratio(target))
            CASE DEFAULT
              STOP 'NO USER TI OPTION YET'
          ENDSELECT
        ENDDO

c...    Loop exit conditions:

      ENDDO


      RETURN
 99   STOP
      END
c
c ====================================================================
c
c...    Set the options specific to the current iteration, which may differ
c       from the option set assigned in the solver initialization routine,
c       particularly on the first solver iteration, or dynamically as a 
c       result of analyzing the present state of the solution:

      SUBROUTINE SetupSolverOptions(count)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none
 
      INTEGER count,ion,target
      LOGICAL, SAVE :: first_bc = .TRUE., set_cnt_bc = .FALSE.

      ion = 1

      IF (log.GT.0) THEN
        WRITE(logfp,*) 
        WRITE(logfp,'(A    )') 'Local control options:'
        WRITE(logfp,'(A,I10)') 'first_bc   = ',first_bc
        WRITE(logfp,'(A,L10)') 'set_cnt_bc = ',set_cnt_bc
      ENDIF

      IF (count.EQ.1) THEN

        machno(0      ,ion) = 1.0D0  ! Store and restore from previous iteration...
        machno(icmax+1,ion) = 1.0D0  ! Should have ION dependence...?

        anl_ic_super(LO) = 0
        anl_ic_super(HI) = icmax + 1
        anl_imaginary    = .FALSE.

        cnt_target       = .TRUE.
        cnt_prescription = .TRUE.
        cnt_particles    = .TRUE.
        cnt_momentum     = .TRUE.
        cnt_energy       = .TRUE.
        cnt_integrate    = .TRUE.

c...    Main loop iteration conidtions:
        cnt_super               = .FALSE.
        cnt_boundary_conditions = .FALSE.

        opt_p_ion_scale = 1      

        DO target = LO, HI

c...      Reassign options if neutrals data not available:
          IF (.NOT.opt%pin_data) THEN

            IF (opt%p_rec(target).EQ.2) opt%p_rec(target) = 1

            IF (opt%p_ion(target).EQ.2) THEN
              opt%p_ion(target) = 3
              opt%p_ion_frac(target) = 0.0
            ENDIF

c...        Sheath limited regime hack (see end of AssignNodeValues):
            IF (opt%p_ion(target).EQ.999) opt%p_ion(target) = 0

            IF (opt%m_mom(target).EQ.2) opt%m_mom(target) = 0
          ELSE
c...        Sheath limited regime hack (see end of AssignNodeValues):
            IF (opt%p_ion(target).EQ.999) opt%p_ion(target) = 2
          ENDIF

c...      When prescribing the ionisation source, scale to sink strength
c         if default source scaling is used (no scaling):
          IF (opt%p_ion     (target).EQ.3.AND.
     .        opt%p_ion_frac(target).EQ.100.0D0) 
     .      opt%p_ion_frac(target) = 0.0

c...      Turn off ionisation source scaling if floating boundary
c         conditions in use:
          IF (opt%bc(target).GE.2) opt_p_ion_scale(target) = 0

          IF (opt%bc(target).LT.0) THEN
            set_cnt_bc = .TRUE.  ! For testing upstream BC's...   *** HACK ***
            opt%bc(target) = -opt%bc(target)
          ENDIF

        ENDDO

      ELSE

        IF (set_cnt_bc.AND.count.EQ.3) cnt_boundary_conditions = .TRUE.  

        IF (cnt_boundary_conditions.AND.
     .      .NOT.cnt_super(LO).AND..NOT.cnt_super(HI)) THEN
          IF (first_bc) THEN
c            first_bc = .FALSE.
            machno(0      ,ion) = 1.0D0
            machno(icmax+1,ion) = 1.0D0
            cnt_target    = .TRUE.
            cnt_momentum  = .TRUE.

            opt%bc          = 3
            opt%p_ion       = 1000
            opt%p_ano       = 1000
            opt%m_mom       = 1000
            opt%m_ano       = 1000
            opt_p_ion_scale = 0

            cnt_energy       = .TRUE.
            cnt_prescription = .TRUE.
            opt%te_ano_psol  = 1000  
            opt%te_ano       = 1000  
            opt%te_ion       = 1000  

            IF (count.EQ.4) THEN
              WRITE(0,*) 'INCREASING TARGET HEAT FLUX'

              CALL IntegrateArray(LO,eneion(1,ion),1,eneint(0,ion))
              WRITE(0,*) '  ENEINT,Qe0:',eneint(icmid,ion),qe(0)

              qe(0) = qe(0) + 0.02 * eneint(icmid,ion)
              eneion(1:icmid,ion) = eneion(1:icmid,ion) * 0.98
              parion(1:icmid,ion) = parion(1:icmid,ion) * 0.98

              CALL IntegrateArray(LO,eneion(1,ion),1,eneint(0,ion))
              WRITE(0,*) '  ENEINT,Qe0:',eneint(icmid,ion),qe(0)
            ENDIF

c            cnt_boundary_conditions = .FALSE.
          ELSE
c            cnt_target    = .TRUE.
c            cnt_momentum  = .TRUE.

c            first_bc = .TRUE.
c            cnt_boundary_conditions = .FALSE.

c            parion(icbnd1(HI):icbnd2(HI),ion) = 
c     .        1.5D0 * parion(icbnd1(HI):icbnd2(HI),ion)  
          ENDIF

      
        ELSE
c...      Turn things off:
          cnt_target       = .FALSE.
          IF (.NOT.set_cnt_bc) cnt_prescription = .FALSE.
          cnt_particles    = .FALSE.
          cnt_momentum     = .FALSE.
          cnt_energy       = .FALSE.
          cnt_integrate    = .FALSE.
        ENDIF 

        IF (cnt_super(LO).OR.cnt_super(HI)) THEN
          cnt_target   = .TRUE.
          cnt_momentum = .TRUE.
        ENDIF

c        cnt_options = .FALSE.
      ENDIF

c...  User control:
      CALL User_SetupSolverOptions(count)

c...  Output:
      IF (log.GT.0) THEN
        WRITE(logfp,*) 
        WRITE(logfp,'(A    )') 'OSM control options:'
        WRITE(logfp,'(A,I10)') 'COUNT        = ',count
        WRITE(logfp,'(A,L10)') 'OPT%PIN_DATA = ',opt%pin_data
        WRITE(logfp,*) 
        WRITE(logfp,'(A,2I10)')   'BC          = ',opt%bc 
        WRITE(logfp,'(A,2I10)')   'SUPER       = ',opt%super
        WRITE(logfp,'(A,2I10)')   'P_REC       = ',opt%p_rec
        WRITE(logfp,'(A,2I10)')   'P_ION       = ',opt%p_ion
        WRITE(logfp,'(A,2I10)')   'P_ION_SCALE = ',opt_p_ion_scale
        WRITE(logfp,'(A,2F10.2)') 'P_ION_FRAC  = ',opt%p_ion_frac
        WRITE(logfp,'(A,2I10)')   'P_ANO       = ',opt%p_ano
        WRITE(logfp,'(A,2I10)')   'M_MOM       = ',opt%m_mom
        WRITE(logfp,'(A,2I10)')   'M_FIT       = ',opt%m_fit
        WRITE(logfp,'(A,2I10)')   'M_ANO       = ',opt%m_ano
        WRITE(logfp,'(A,2I10)')   'M_ANO_DIST  = ',opt%m_ano_dist
        WRITE(logfp,'(A,2F10.2)') 'M_ANO_EXP   = ',opt%m_ano_exp
        WRITE(logfp,'(A,2I10)')   'TE_REC      = ',opt%te_rec
        WRITE(logfp,'(A,2I10)')   'TE_ION      = ',opt%te_ion
        WRITE(logfp,'(A,2I10)')   'TE_ANO      = ',opt%te_ano
        WRITE(logfp,'(A,2I10)')   'TE_ANO_PSOL = ',opt%te_ano_psol
        WRITE(logfp,'(A,2I10)')   'TE_CONV     = ',opt%te_conv
        WRITE(logfp,'(A,2I10)')   'TI          = ',opt%ti
        WRITE(logfp,'(A,2I10)')   'TI_REC      = ',opt%ti_rec
        WRITE(logfp,'(A,2I10)')   'TI_ION      = ',opt%ti_ion
        WRITE(logfp,'(A,2I10)')   'TI_ANO      = ',opt%ti_ano
        WRITE(logfp,'(A,2I10)')   'TI_ANO_PSOL = ',opt%ti_ano_psol
        WRITE(logfp,'(A,2I10)')   'TI_CONV     = ',opt%ti_conv
        WRITE(logfp,'(A,2F10.2)') 'TI_RATIO    = ',opt%ti_ratio
        WRITE(logfp,'(A,2I10)')   'TI_EQUIL    = ',opt%ti_equil
        WRITE(logfp,*) 
      ENDIF

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE AssignTargetConditions
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      REAL*8 GetNodePressure

      INTEGER ion,target,ic1,ic2,ic,in
      REAL*8  net(3),integral(10),cs,source(icmax),
     .        srcint(0:icmax),fact

      LOGICAL, SAVE :: firstcall = .TRUE.  ! TEMP

      DO ion = 1, nion
        IF (iontype(ion).NE.ITY_FLUID) CYCLE
        net = 0.0D0
        DO target = LO, HI
           ic = ictarg(target)
           in = intarg(target)

c...      Set target values:
c         --------------------------------------------------------------
          SELECTCASE (opt%bc(target))
            CASE (1)
c...          Target conditions specified:

              te(ic)     = DBLE(node(in)%te)  ! Add some checks...
              ti(ic,ion) = te(ic) * DBLE(opt%ti_ratio(target))  ! Need to check options: Ti from input file or ratio, etc...

              IF     (node(in)%jsat(ion).NE.0.0) THEN
                fact = area(ic) / ECH
                isat(ic,ion) = -DBLE(node(in)%jsat(ion)) * fact

c        WRITE(logfp,*) 'ISAT0:',isat(ic,ion),node(in)%jsat(ion),fact
c                isat(ic,ion) = -DBLE(ABS(node(in)%jsat(ion))) * fact
              ELSEIF (node(in)%ne       .NE.0.0) THEN         
                WRITE(0,*) ' NODE IN=',in
                STOP 'NE TARGET FLUX NOT READY'
              ELSEIF (node(in)%pe       .NE.0.0) THEN        
                STOP 'PE TARGET FLUX NOT READY'
              ELSE
                CALL ER('AssignTargetConditions','Target particle '//
     .                  'flux not specified',*99)
              ENDIF

              cs = DSQRT( (te(ic) + ti(ic,ion)) * ECH / mi(ion) )  ! Needs improvement... CalcCs

c             IPP/08 Krieger - SUN compiler insists on parentheses around -1
              vi(ic,ion) = cs * machno(ic,ion) * tsign(target) *
     .                     (-1.0D0) * DSIGN(1.0D0,isat(ic,ion))

c              vi(ic,ion) = cs * machno(ic,ion) * tsign(target) *
c     .                     -1.0D0 * DSIGN(1.0D0,isat(ic,ion))

              ne(ic)     = DABS(isat(ic,ion) / vi(ic,ion))
              ni(ic,ion) = ne(ic)
              pe(ic)     = ne(ic) * te(ic) * ECH 
              pi(ic,ion) = ni(ic,ion) * 
     .                     (ti(ic,ion) * ECH + mi(ion) * vi(ic,ion)**2)

            CASE (2)
c...          Target particle and momentum from upstream cross-field 
c             and volume sources:

              te(ic)     = DBLE(node(in)%te)  ! Add some checks...
              ti(ic,ion) = te(ic) * DBLE(opt%ti_ratio(target))  ! Not good, need to check options

              cs = DSQRT( (te(ic) + ti(ic,ion)) * ECH / mi(ion) )  ! Needs improvement... CalcCs
              vi(ic,ion) = cs * machno(ic,ion) * tsign(target)

              isat(ic,ion) =  0.0D0
              ne(ic)       = -1.0D0
              ni(ic,ion)   = -1.0D0
              pe(ic)       = -1.0D0
              pi(ic,ion)   = -1.0D0

            CASE (3)
c...          Floating jsat and Te at the target:

c...TEMP: ...need to get this from elsewhere...
              IF (firstcall) THEN
                qe(ic) = 3.0D0 * (ECH*te(ic)) * isat(ic,ion)

                WRITE(0,*) 'QE:',qe(ic),te(ic),isat(ic,ion)*ECH
                IF (target.EQ.HI) qe(ic) = -1.0D0 * qe(ic)
                IF (target.EQ.HI) firstcall = .FALSE.
              ENDIF
c              STOP 'sdgsdgsdg'

              isat(ic,ion) =  0.0D0
              ne(ic)       = -1.0D0
              ni(ic,ion)   = -1.0D0
              vi(ic,ion)   = -1.0D0
              pe(ic)       = -1.0D0
              pi(ic,ion)   = -1.0D0
              te(ic)       = -1.0D0
              ti(ic,ion)   = -1.0D0

            CASEDEFAULT
              CALL ER('AssignTargetConditions','Unknown option',*99)
          ENDSELECT

        ENDDO  ! End of target loop
      ENDDO  ! End of ION loop

      WRITE(logfp,*) 'ISAT1:',ion,isat(ictarg(LO),1),
     .                            isat(ictarg(HI),1)

      RETURN
 99   WRITE(0,*) '  TARGET=',target
      WRITE(0,*) '  OPT%BC=',opt%bc(target)
      STOP
      END
c
c ====================================================================
c
c... 
c
      SUBROUTINE EvaluateFluidSolution
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none
 
      INTEGER target,ictarget(2),ion
      REAL*8  delta(2),adjust(2)

      DATA delta, adjust /0.5D0, 0.5D0, 0.0D0, 0.0D0/

      SAVE

      ion = 1

      ictarget(LO) = 0
      ictarget(HI) = icmax + 1

      DO target = LO, HI

c...    Check for near-target sonic transition:
        SELECTCASE (opt%super(target))
          CASE(0)
c           Do nothing:
          CASE(1)
c           Delta target Mach number and iterate:
            IF (anl_imaginary(target)) THEN
              IF     (adjust(target).EQ.0.0D0) THEN
                adjust(target) =  0.1D0
              ELSEIF (adjust(target).LT.0.0D0) THEN
                adjust(target) = -0.5D0 * adjust(target)
              ENDIF
              cnt_super(target) = .TRUE. 
            ELSE
              IF     (adjust(target).GT.0.0D0) THEN
                adjust(target) = -0.5D0 * adjust(target)
                cnt_super(target) = .TRUE. 
              ELSEIF (DABS(adjust(target)).LT.1.0D-3) THEN
c...            Done, reset:
                delta(target) = 0.5D0
                adjust(target) = 0.0D0
                cnt_super(target) = .FALSE. 
              ENDIF
            ENDIF

            IF (cnt_super(target)) THEN
              delta(target) = delta(target) * (1.0D0 + adjust(target))
              machno(ictarget(target),ion) = 1.0D0 + delta(target)
            ENDIF

            IF (logfp.GT.0)     
     .        WRITE(logfp,'(A,I4,L2,I4,3F12.6,L2)') 'SUPER: ',
     .          target,anl_imaginary(target),
     .          anl_ic_super(target),delta(target),adjust(target),
     .          machno(ictarget(target),ion),cnt_super(target)

          CASEDEFAULT
            STOP 'NO USER ROUTINES YET'
        ENDSELECT



      ENDDO

      RETURN
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE SOL28_V4(tube1,icmax1,cell1,pin1,fluid1,
     .                    ref_tube1,ref_nion1,ref_icmax1,ref_fluid1,  ! Do these need the "1" business?
     .                    nnode1,mnode1,node1,nion1,opt_global)
      USE mod_sol28_params
      USE mod_sol28
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER, INTENT(IN) :: icmax1,nnode1,mnode1,nion1,
     .                       ref_nion1,ref_icmax1
      TYPE(type_tube   ) :: tube1
      TYPE(type_cell   ) :: cell1     (icmax1)
      TYPE(type_neutral) :: pin1      (icmax1,nion1)
      TYPE(type_fluid  ) :: fluid1    (icmax1,nion1)
      TYPE(type_tube   ) :: ref_tube1
      TYPE(type_fluid  ) :: ref_fluid1(ref_icmax1,ref_nion1)
      TYPE(type_node   ) :: node1     (nnode1)
      TYPE(type_options_osm ), INTENT(IN) :: opt_global

      INTEGER count,ion,ic
      LOGICAL cont

c...  Halt the code since data is almost certainly not being
c     passed properly, needs fixing but likely this is not an
c     issue for some time:
      IF (nion1.NE.1) CALL ER('SOL28_V4','NION not equal to 1',*99)

c...  Data transfer:  ...don't really like this... other way of doing it?
      nion = nion1
      tube = tube1
      icmax = icmax1  
      cell(1:icmax) = cell1(1:icmax)
      DO ion = 1, nion
        pin  (1:icmax,ion) = pin1  (1:icmax,ion)
        fluid(1:icmax,ion) = fluid1(1:icmax,ion)
      ENDDO
      IF (ref_icmax1.GT.1) THEN
        ref_tube = ref_tube1
        IF (icmax1.NE.ref_icmax1) 
     .    CALL ER('SOL28_V4','Reference data does not match tube',*99)
        ref_nion = ref_nion1
        WRITE(0,*) 'ref_nion,ref_icmax:',
     .    ref_nion,ref_icmax1
        DO ion = 1, ref_nion
          ref_fluid(1:icmax,ion) = ref_fluid1(1:icmax,ion)
        ENDDO
      ELSE
        ref_nion = 0
      ENDIF

      
      nnode = nnode1
      mnode = mnode1
      node(1:nnode) = node1(1:nnode) 

      opt = opt_global
      log = opt%log
      logfp = opt%logfp

c      WRITE(0,*) '-->',opt%bc(1:2)

c...  (Some assignments):
      iontype(1) = ITY_FLUID

      icmid = node(mnode)%icell

      icbnd1(LO) = 1
      icbnd2(LO) = icmid
      icbnd1(HI) = icmid + 1
      icbnd2(HI) = icmax

c      WRITE(0,*) 'ICBND:',icbnd1
c      WRITE(0,*) 'ICBND:',icbnd2

      ictarg(LO) = 0
      ictarg(HI) = icmax + 1
      intarg(LO) = 1
      intarg(HI) = nnode

      tsign(LO) = -1.0D0
      tsign(HI) =  1.0D0

      sfor(0      ) = 0.0D0
      sbak(0      ) = tube%smax
      sfor(1:icmax) = DBLE(            cell(1:icmax)%s)
      sbak(1:icmax) = DBLE(tube%smax - cell(1:icmax)%s)
      sfor(icmax+1) = tube%smax
      sbak(icmax+1) = 0.0D0
      smax = DBLE(tube%smax)

      sdelta(1:icmax) =DBLE(cell(1:icmax)%sbnd(2)-cell(1:icmax)%sbnd(1))
      area(0:icmax+1) = 1.0D0
      vol (1:icmax  ) = DBLE(cell(1:icmax)%ds)

      ai(1) = 2.0D0
      mi(1) = ai(1) * AMU
      zi(1) = 1.0D0
      ci(1) = 1.0D0

c...  Basic initialization:
      cont = .TRUE.
      count = 0

c...  Control:
      cnt_options = .TRUE.

      IF (log.GT.0) THEN
        WRITE(logfp,*) 'M,NNODE:',mnode,nnode
        WRITE(logfp,*) '  1    :',node(1    )%s ,node(1    )%te,
     .                            node(1    )%ne,node(mnode)%pe
        WRITE(logfp,*) '  mnode:',node(mnode)%s ,node(mnode)%te,
     .                            node(mnode)%ne,node(mnode)%pe
        WRITE(logfp,*) '  nnode:',node(nnode)%s ,node(nnode)%te,
     .                            node(nnode)%ne,node(nnode)%pe
      ENDIF


      ne = 0.0D0
      pe = 0.0D0
      te = 0.0D0

      ni = 0.0D0
      pi = 0.0D0
      vi = 0.0D0
      ti = 0.0D0

      momvol = 0.0D0

      qcond = 0.0D0
      qconv = 0.0D0

c...  Main loop:
      DO WHILE (cont)
c...    Loop counter:
        count = count + 1

        IF (log.GT.0) WRITE(logfp,*) '    COUNT:',tube%ir,count

c...    Set options for solver iteration:
c       ----------------------------------------------------------------
        IF (cnt_options) THEN
          CALL SetupSolverOptions(count)
        ENDIF

c...    Compile/modify ion quantities from kinetic codes:
c       ----------------------------------------------------------------
        IF (.TRUE.) THEN
          CALL ProcessKineticIons  ! Only needs to be done once?
        ENDIF

c...    Boundary conditions:
c       ----------------------------------------------------------------
        IF (cnt_target) THEN 
          CALL AssignTargetConditions
        ENDIF

  
        WRITE(logfp,*) 'ISAT2:',ion,isat(ictarg(LO),1),
     .                              isat(ictarg(HI),1)

c...    Set prescribed quantities:
c       ----------------------------------------------------------------
        IF (cnt_energy) THEN
          CALL AssignEnergySources   ! This goes below... for floating boundary conditions... Te,i assigned above
        ENDIF

        IF (cnt_prescription) THEN
c          te = 10.0D0
          CALL PrescribeTemperatureProfiles  ! This goes in the solve fluid equations routine I think...
        ENDIF

c...    Process source terms:
c       ----------------------------------------------------------------
        IF (cnt_particles) THEN
          CALL AssignParticleSources
        ENDIF

        IF (cnt_momentum) THEN
          CALL AssignMomentumSources
        ENDIF

        IF (cnt_particles.OR.cnt_momentum) THEN
          CALL IntegrateSources(2)
        ENDIF

c...    In the case of floating boundary conditions:
c       ----------------------------------------------------------------
        IF (cnt_boundary_conditions) THEN  ! Change to cnt_floating_target...
c. LEFT OFF: target jsat and Te are stable
c     -need improved QE calculation maintenance, and need to update the value when
c      the volume source on the half-ring changes... simulate artificially first
c     -need logic in UpdateTemperatureProfiles to evolve upstream parameters to give correct
c      target values...
c     -need an exit condition for this?  or just during testing?


          CALL UpdateTemperatureProfiles 
        ENDIF

c...    Solve fluid equations:
c       ----------------------------------------------------------------
        CALL SolveFluidEquations

c...    Evaluate solution:
c       ----------------------------------------------------------------
        CALL EvaluateFluidSolution


c...    Main iteration loop exit conditions:
c       ----------------------------------------------------------------
        IF     (count.LE.3) THEN
          cont = .TRUE.
        ELSEIF (count.GT.50) THEN
          WRITE(logfp,*) 'SORRY, GIVING UP... NO CONVERGENCE'
          cont = .FALSE.
        ELSEIF ((cnt_super(LO).OR.cnt_super(HI))) THEN
c        ELSEIF (.FALSE..AND.(cnt_super(LO).OR.cnt_super(HI))) THEN
c...      Near-target sonic transition being processed:
          cont = .TRUE. 
        ELSEIF (.FALSE.) THEN
c...      Source/flux adjustment required: 
          cont = .TRUE. 
        ELSEIF (.FALSE.) THEN
c...      Ensure sufficient iterations for converged solution: 
          cont = .TRUE. 
        ELSEIF (cnt_boundary_conditions) THEN
c...      Establish cross-field source terms for upstream boundary
c         conditions:
          cont = .TRUE. 
c        ELSEIF (count.EQ.1) THEN   
c          cont = .TRUE.
        ELSE
          cont = .FALSE.
        ENDIF

      ENDDO ! End of main loop


c...  Assign plasma quantities to tube arrays:
      ion = 1

      ic = ictarg(LO)
      WRITE(logfp,*) 'JSAT5:',tube%jsat(LO,ion),ech
      tube%jsat  (LO,ion) = -SNGL(isat  (ic,ion)) * ECH  ! Change units in the solver to recover A m-2... what is MKS?
      WRITE(logfp,*) 'JSAT5:',tube%jsat(LO,ion)
      tube%ne    (LO)     = SNGL(ne    (ic))
      tube%te    (LO)     = SNGL(te    (ic))
      tube%ni    (LO,ion) = SNGL(ni    (ic,ion))
      tube%vi    (LO,ion) = SNGL(vi    (ic,ion))
      tube%machno(LO)     = SNGL(machno(ic,ion))  ! Should this have an ION index or not...
      tube%ti    (LO,ion) = SNGL(ti    (ic,ion))

      ic = ictarg(HI)
      WRITE(logfp,*) 'JSAT5:',tube%jsat(HI,ion),ech
      tube%jsat  (HI,ion) = -SNGL(isat  (ic,ion)) * ECH
      WRITE(logfp,*) 'JSAT5:',tube%jsat(HI,ion),ech
      tube%ne    (HI)     = SNGL(ne    (ic))
      tube%te    (HI)     = SNGL(te    (ic))
      tube%ni    (HI,ion) = SNGL(ni    (ic,ion))
      tube%vi    (HI,ion) = SNGL(vi    (ic,ion))
      tube%machno(HI)     = SNGL(machno(ic,ion))  ! Should this have an ION index or not...
      tube%ti    (HI,ion) = SNGL(ti    (ic,ion))

      DO ion = 1, nion
        fluid(1:icmax,ion)%ne = SNGL(ne(1:icmax))
        fluid(1:icmax,ion)%te = SNGL(te(1:icmax))
        fluid(1:icmax,ion)%ni = SNGL(ni(1:icmax,ion))
        fluid(1:icmax,ion)%vi = SNGL(vi(1:icmax,ion))
        fluid(1:icmax,ion)%ti = SNGL(ti(1:icmax,ion))

        fluid(1:icmax,ion)%parion =  SNGL(parion(1:icmax,ion))
        fluid(1:icmax,ion)%parrec = -SNGL(parrec(1:icmax,ion))
        fluid(1:icmax,ion)%parano =  SNGL(parano(1:icmax,ion))

c        WRITE(0,*) 'STRANGE:',fluid(1:icmax,ion)%parion

        fluid(1:icmax,ion)%momvol =  SNGL(momvol(1:icmax,ion))
        fluid(1:icmax,ion)%momano =  SNGL(momano(1:icmax,ion))

        fluid(1:icmax,ion)%enerec =  SNGL(enerec(1:icmax,ion))
        fluid(1:icmax,ion)%eneion =  SNGL(eneion(1:icmax,ion))

        fluid(1:icmax,ion)%eneano =  SNGL(eneano(1:icmax    ))  ! Some inconsistency in the use of "ene"...

        fluid(1:icmax,ion)%parusr = SNGL(parusr(1:icmax,ion))
        fluid(1:icmax,ion)%momusr = SNGL(momusr(1:icmax,ion))
        fluid(1:icmax,ion)%eneusr = SNGL(eneusr(1:icmax))
        fluid(1:icmax,ion)%eniusr = SNGL(eniusr(1:icmax,ion))

        fluid(1:icmax,ion)%parsrc =  SNGL(parsrc(1:icmax,ion))
        fluid(1:icmax,ion)%momsrc =  SNGL(momsrc(1:icmax,ion))
        fluid(1:icmax,ion)%enesrc =  SNGL(enesrc(1:icmax,1))
      ENDDO

c...  Pass data back through global arrays (ugly):
      tube1 = tube
      DO ion = 1, nion
        fluid1(1:icmax,ion) = fluid(1:icmax,ion)
      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE MainLoop(itstart,itend,ikopt)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      INTEGER itstart,itend,ikopt,count

      INTEGER itube,cind1,cind2,ref_cind1,ref_cind2,ref_itube
      LOGICAL cont

      INTEGER nnode,mnode
      TYPE(type_node) :: node(S28_MAXNNODE)
      TYPE(type_options_osm) :: opt_tube

      IF (log.GT.0) THEN
        WRITE(logfp,*)
        WRITE(logfp,*) 'RANGE:',itstart,itend,grid%isep
        WRITE(logfp,*) 'TE   :',tube(itstart)%te(LO),
     .                          tube(itstart)%te(HI)
      ENDIF

      opt%cosm = 0

      count = 0
      cont  = .TRUE.
      DO WHILE (cont) 
        count = count + 1
        cont  = .FALSE.
 

c...    Calcualte derived quantities (efield, drifts, etc.):
c        CALL CalculateEfield
c        CALL CalculateGradients
c        CALL CalculateDrifts

        DO itube = itstart, itend

c...      Setup solver options:
c         CALL SetupLocalOptions(itube,opt)

          opt_tube = opt

          IF (log.GT.0) THEN
            WRITE(logfp,'(A   )') '------------------------------------'
            WRITE(logfp,'(A,I6)') ' SOLVING TUBE: ',itube
            WRITE(logfp,'(A,I6)') ' CFLUKIN     : ',opt_tube%cflukin
            WRITE(logfp,'(A,I6)') ' COSM        : ',opt_tube%cosm
            WRITE(logfp,'(A   )') '------------------------------------'
          ENDIF
          
c...      Assign solution parameter nodes:
          IF (opt%s28mode.EQ.4.1) THEN 
            CALL AssignNodeValues_New(itube,nnode,mnode,node,opt_tube)
          ELSE
            CALL AssignNodeValues_Legacy(itube,nnode,mnode,node)
          ENDIF

c...      Identify cell based data to transfer to solver:
          cind1 = tube(itube)%cell_index(LO)
          cind2 = tube(itube)%cell_index(HI)

          ref_itube = MIN(ref_ntube,itube)
          ref_cind1 = MIN(ref_nfluid,cind1)
          ref_cind2 = MIN(ref_nfluid,cind2)

 
c...      Calculate plasma solution:
          IF (.FALSE.) THEN
c...        MPI:
          ELSE
            IF (tube(itube)%type.EQ.GRD_CORE) THEN
              CALL AssignCorePlasma(itube,nnode,mnode,node)
            ELSE

              CALL SOL28_V4(tube(itube),tube(itube)%n,                   ! Likely need to pass index range for 
     .                      cell (cind1:cind2),                          ! each array individually in case
     .                      pin  (cind1:cind2,nion),                     ! some have only nominal allocations, i.e.
     .                      fluid(cind1:cind2,nion),                     ! they are not in use...
     .                      ref_tube(ref_itube),
     .                      ref_nion,ref_cind2-ref_cind1+1,
     .                      ref_fluid(ref_cind1:ref_cind2,ref_nion),
     .                      nnode,mnode,node,nion,opt_tube)                   ! Also: pass local options

c              WRITE(0,*) 'TUBE:PSOL',tube(itube)%eneano(1:2)

            ENDIF
          ENDIF 

        ENDDO  ! Tube loop

c...    Conditions for iteration:
        CALL User_MainLoop(cont,count)

      ENDDO ! Iteration loop

c      CALL OutputData(logfp,'FLUID SOLUTION COMPLETE')

      RETURN
 99   STOP
      END

