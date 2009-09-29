c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE AssignParticleSources
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER ion, target, ic1, ic2
      REAL*8  source(icmax),decay

      cnt_integrate = .TRUE.
 
      IF (log.GT.1) WRITE(logfp,*) 'ASSIGNING PARTICLE SOURCES'

      DO ion = 1, nion
        IF (iontype(ion).NE.ITY_FLUID) CYCLE
        DO target = LO, HI
          ic1 = icbnd1(target)
          ic2 = icbnd2(target)

c           WRITE(0,*) 'IC1,2:',ic1,ic2

c...      Volume:
c         ----------------------------------------------------------------

c...      Recombination:
          source = 0.0D0
          SELECTCASE (opt%p_rec(target))
            CASE (1000)
c             Preserve from previous iteration:               
              source(ic1:ic2) = parrec(ic1:ic2,ion)
            CASE (0)
c             None:          
            CASE (1)
c             Reference plasma:
              source(ic1:ic2) = -DBLE(ref_fluid(ic1:ic2,ion)%parrec) 
            CASE (2)
c             PIN:          
              source(ic1:ic2) = -1.0D0 * DBLE(pin(ic1:ic2,ion)%rec) 
c            CASE (3)
c              Calculated...:
            CASEDEFAULT
              CALL User_VolumeParRecSource(target,source)
          ENDSELECT         
          parrec(ic1:ic2,ion) = source(ic1:ic2)

c...      Ionisation:
          source = 0.0D0
          SELECTCASE (opt%p_ion(target))
            CASE (1000)
c             Preserve (from a previous iteration):
              source(ic1:ic2) = parion(ic1:ic2,ion)               
            CASE (0)
c             None:          
            CASE (1)
c             Reference plasma:
              source(ic1:ic2) = DBLE(ref_fluid(ic1:ic2,ion)%parion) 
c              WRITE(logfp,*) 'parion:',source(ic1:ic2)
            CASE (2)
c             PIN:
              source(ic1:ic2) = DBLE(pin(ic1:ic2,ion)%ion) 
              WRITE(logfp,*) 'parion:',source(ic1:ic2)
            CASE (3)
c             Prescribed:
c              IF (target.EQ.LO) THEN 
c                CALL SpecifyDistribution(target,-1,0,2,10.D0,source)  ! PROMOTES SUPERSONIC TARGETS 
c              WRITE(logfp,*) '--> IONSRC ',target
c                CALL SpecifyDistribution(target,-2,0,2,DBLE(opt%p_ion_frac(target)),source) 
c              ELSE
              decay = DBLE(opt%p_ion_exp(target))
              CALL SpecifyDistribution(target,-2,0,2,decay,source)
c              ENDIF
c              WRITE(logfp,*) '--> DONE'
            CASEDEFAULT                                            
              CALL User_VolumeParIonSource(target,source)          
          ENDSELECT         
          parion(ic1:ic2,ion) = source(ic1:ic2)


c...      Cross-field:
c         ----------------------------------------------------------------


c...      Anomalous:
c         ----------------------------------------------------------------
          source(ic1:ic2) = 0.0D0
          SELECTCASE (opt%p_ano(target))
            CASE(1000)
c...          Preserve (from a previous iteration):
              source(ic1:ic2) = parano(ic1:ic2,ion)               
            CASE (0)
c...          Do nothing:
            CASE (1)
c             Reference plasma:
              source(ic1:ic2) = DBLE(ref_fluid(ic1:ic2,ion)%parano)
            CASE (2)
c...          Half ring:
              CALL SpecifyDistribution(target,0,0,1,0.0D0,source)
            CASE (3)
c...          Full ring:
              IF (opt%bc(LO).NE.1) 
     .          CALL ER('AssignParticleSources','Full ring '//    ! Move this check to SetupSolverOptions...
     .                  'anomalous cross-field flux '//
     .                  'specification requires conditions at '// 
     .                  'both targets to be specified',*99)       

              CALL SpecifyDistribution(FULL,0,0,1,0.0D0,source)   ! This will be moved to the ConserveParticles
            CASEDEFAULT                                           ! routine and operate as node interpolation does
              STOP 'USER ROUTINE NOT READY: P_ANO'                ! for momentum?  But this doesn't seem rights... (pain)
              CALL User_VolumeParAnoSource(target,source)         ! Need to add P_ANO_DIST, etc. as for M_ANO
          ENDSELECT
          parano(ic1:ic2,ion) = source(ic1:ic2)

c...      Generic user particle source specification:
c         ----------------------------------------------------------------
          source = 0.0D0
          CALL User_ParticleSource(ion,target,source)
          parion(ic1:ic2,ion) = parion(ic1:ic2,ion) + source(ic1:ic2)  ! *HACK* source needs to go somewhere else...

        ENDDO
      ENDDO

c...
      CALL ConserveParticles

c...
      CALL IntegrateSources(1)

      IF (log.GT.1) WRITE(logfp,*) 'ASSIGNING PARTICLE SOURCES: DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE ConserveParticles
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none


      INTEGER ion,target,ict,ic1,ic2 , i1 ! TEMP
      REAL*8  net(3),integral(10),fact,intrec,diff,frac,scale


      DO ion = 1, nion
        IF (iontype(ion).NE.ITY_FLUID) CYCLE
        net = 0.0D0
        DO target = LO, HI
          ict = ictarg(target)
          ic1 = icbnd1(target)
          ic2 = icbnd2(target)

c...      Scale volume recombination source?
c         --------------------------------------------------------------

c...      Scale volume ionisation source:
c         --------------------------------------------------------------

          SELECTCASE (opt_p_ion_scale(target))
            CASE (0)
c             Do nothing:
            CASE (1)
c             Enforce minimum particle balance on each half ring. Get
c...          particle sink for the half ring and limit the ionisation
c             source to some fraction of the total sink:

              scale = 1.0D0

              frac = DBLE(opt%p_ion_frac(target))

              CALL IntegrateArray(target,parrec(1,ion),0,integral(1))
              CALL IntegrateArray(target,parion(1,ion),0,integral(2))
              net(target) = -1.0D0 * (isat(ict,ion) + integral(1))

              diff = (integral(2) - net(target)) / net(target)
              IF ((diff.GT. 1.0D0 * frac).OR.
     .            (diff.LT.-1.0D0 * frac)) THEN 
                scale = (net(target) * (1.0D0 - frac)) / integral(2)
              ENDIF

              parion(ic1:ic2,ion) = parion(ic1:ic2,ion) * scale
c...TEMP: for testing...
              IF (.FALSE..AND.opt%p_ano(target).EQ.3) THEN
                parion(ic1:ic2,ion) = parion(ic1:ic2,ion) * 0.9
                WRITE(0,*) 'WARNING: SCALING PARION SOURCE FOR TESTING'
              ENDIF

            CASEDEFAULT
              CALL ER('ConserveParticles','Unknown P_ION option',*99)
          ENDSELECT

c...      Scale anomalous source:
c         --------------------------------------------------------------
          SELECTCASE (opt%p_ano(target))
            CASE (1000)
            CASE (0)
c...          Do nothing (set to 0.0 in AssignParticleSources):
            CASE (1)
c             Don't overwrite what was set in AssignParticleSources from a 
c             previous plasma solution:
            CASE (2)
c             Half ring:
              CALL IntegrateArray(target,parrec(1,ion),0,integral(1))
              CALL IntegrateArray(target,parion(1,ion),0,integral(2))
              net(target) = isat(ict,ion) + integral(1) + integral(2)
              parano(ic1:ic2,ion) = -parano(ic1:ic2,ion) * net(target)
            CASE (3)
c             Full ring:

              IF (opt%bc(LO).NE.1) 
     .          CALL ER('AssignParticleSources','Full ring '//    ! Move this check to SetupSolverOptions...
     .                  'anomalous cross-field flux '//
     .                  'specification requires conditions at '// 
     .                  'both targets to be specified',*99)

              IF (target.EQ.HI) THEN
                CALL IntegrateArray(FULL,parrec(1,ion),0,integral(1))
                CALL IntegrateArray(FULL,parion(1,ion),0,integral(2))
                net(FULL) = isat(ictarg(LO),ion) + 
     .                      isat(ictarg(HI),ion) + 
     .                      integral(1) + integral(2)
                parano(1:icmax,ion) = -parano(1:icmax,ion) * net(FULL)
              ENDIF
            CASEDEFAULT
              CALL ER('ConserveParticles','Unknown P_ANO option',*99)
          ENDSELECT

        ENDDO  ! End of TARGET loop
      ENDDO  ! End of ION loop

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE AssignMomentumSources
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER ion, target, ic1, ic2
      REAL*8  source(icmax)

      IF (log.GT.1) WRITE(logfp,*) 'ASSIGNING MOMENTUM SOURCES'
  
c...  Set everything to zero?  No, because for performance optimization I really
c     need to keep the number of non-necessary re-assignments to a minimum, 
c     so need to blank everything to start and then assign on the 1st iteration only...

      cnt_integrate = .TRUE.

      DO ion = 1, nion
        IF (iontype(ion).NE.ITY_FLUID) CYCLE
        DO target = LO, HI
          ic1 = icbnd1(target)
          ic2 = icbnd2(target)

c...      Volume:
c         ----------------------------------------------------------------

c...      Neutral friction:
          source = 0.0D0
          SELECTCASE (opt%m_mom(target))
            CASE (1000)
c             Preserve:
              source(ic1:ic2) = momvol(ic1:ic2,ion)
            CASE (0)
c             None:          
              source(ic1:ic2) = 0.0D0
            CASE (1)
c             Reference plasma:          
              source(ic1:ic2) = DBLE(ref_fluid(ic1:ic2,ion)%momvol) 
            CASE (2)
c             PIN:          
              source(ic1:ic2) = DBLE(pin(ic1:ic2,ion)%mom)  ! add scaling factor..?
c            CASE (3)
c             CX based estimate (PIN data required):
            CASEDEFAULT
              CALL ER('MomentumConservation','Unrecognized option',*99)
          ENDSELECT         
          momvol(ic1:ic2,ion) = source(ic1:ic2)

c...      Cross-field:
c         ----------------------------------------------------------------

c...      User:
c         ----------------------------------------------------------------
          source = 0.0D0
          CALL User_MomentumSource(ion,target,source)
          momvol(ic1:ic2,ion) = momvol(ic1:ic2,ion) + source(ic1:ic2)

        ENDDO
      ENDDO

c...
      CALL ConserveMomentum

c...
      CALL IntegrateSources(2)

      IF (log.GT.1) WRITE(logfp,*) 'ASSIGNING MOMENTUM SOURCES: DONE'

      RETURN
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE ConserveMomentum
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      REAL*8 GetNodePressure

      INTEGER ion,target,ic1,ic2,ic,node1,node2,i1,i2,snode,enode,
     .        inode
      REAL*8  net,integral(10),p1,p2,source(icmax),pmatch,
     .        srcint(0:icmax)

      DO ion = 1, nion
        IF (iontype(ion).NE.ITY_FLUID) CYCLE
c...    Anomalous pressure source to enforce pressure evolution between
c       interpolation nodes:
c       ----------------------------------------------------------------
        node1 = 1
        DO WHILE(node1.LT.nnode)

          IF (node1.GT.1.AND.
     .        node(node1)%ne.LE.0.0.AND.
     .        node(node1)%pe.LE.0.0) THEN
            node1 = node1 + 1
            CYCLE
          ENDIF

          DO node2 = node1+1, nnode
            IF (node2.EQ.nnode.OR.
     .          node(node2)%ne.GT.0.0.OR.
     .          node(node2)%pe.GT.0.0) EXIT                 
          ENDDO

          IF (node1.LT.mnode) THEN
            target = LO
          ELSE
            target = HI
          ENDIF

          ic1 = MAX(node(node1)%icell,1)      ! icbnd1(target)
          ic2 = MIN(node(node2)%icell,icmax)  ! icbnd2(target)

c          WRITE(0,*) 'NODES:',node1,node2,ic1,ic2,target

          DO inode = node1+1, node2-1
            IF (node(inode)%ne.EQ.-1.0.OR.
     .          node(inode)%pe.EQ.-1.0) THEN
              IF (inode.LE.mnode) ic1 = node(inode)%icell
              IF (inode.GE.mnode) ic2 = node(inode)%icell
            ENDIF
          ENDDO
          IF (node1.GT.1) ic1 = ic1 + 1  ! Avoids double counting of common node

          source = 0.0D0
          SELECTCASE (opt%m_ano(target))
            CASE (1000)
c             Preserve:
              source(ic1:ic2) = momano(ic1:ic2,ion)
            CASE (0)
c             None:          
            CASE (1)
c...          Reference plasma:
              source(ic1:ic2) = DBLE(ref_fluid(ic1:ic2,ion)%momano)
            CASE (2)
c...          Check for node density/pressure specifications:

c              WRITE(0,*) 'NODES:',node1,node2,ic1,ic2,target

              p1 = GetNodePressure(node1,ion)                 ! Should I add a 1/2 cell momentum
              p2 = GetNodePressure(node2,ion)                 ! array, so that densities are always

              CALL SpecifyDistribution(target,ic1,ic2,
     .                                 opt%m_ano_dist(target),
     .                            DBLE(opt%m_ano_exp (target)),source)


c              WRITE(0,*) 'DIST:',source(ic1:ic2)

              CALL IntegrateArray(FULL,momvol(1,ion),1,srcint)   ! Not debugged as yet...

              IF     (node1.EQ.1.AND.node2.EQ.nnode) THEN
                srcint(0) = srcint(TOTAL)
              ELSEIF (node1.EQ.1) THEN
                srcint(0) = srcint(ic2)
              ELSEIF (node2.EQ.nnode) THEN
                srcint(0) = srcint(TOTAL) - srcint(ic1)
              ELSE                
                srcint(0) = srcint(ic2) - srcint(ic1)
              ENDIF

              source(ic1:ic2) = source(ic1:ic2) * (p2 - p1 - srcint(0))  ! TOTAL=0

c              WRITE(0,*) 'SRC:',p1,p2,srcint(0)
c              WRITE(0,*) source(ic1:ic2)

            CASE DEFAULT
              STOP 'NO USER MOM READY'
              CALL User_MomentumSource(target,source)
          ENDSELECT         
          momano(ic1:ic2,ion) = source(ic1:ic2)

          node1 = node2

        ENDDO  ! End of node loop
      ENDDO  ! End of ION loop

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE AssignEnergySources
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER ion, target, ic1, ic2
      REAL*8  source(icmax)

      IF (log.GT.1) WRITE(logfp,*) 'ASSIGNING ENERGY SOURCES'

      DO ion = 1, nion
        IF (iontype(ion).NE.ITY_FLUID) CYCLE
        DO target = LO, HI
          ic1 = icbnd1(target)
          ic2 = icbnd2(target)

c...      Volume:
c         ----------------------------------------------------------------

c...      PINQe:
          source = 0.0D0
          SELECTCASE (opt%te_rec(target))
c            CASE (1000)
c             Preserve from previous iteration:               
c              source(ic1:ic2) =  ! parrec(ic1:ic2,ion)
c            CASE (0)
c             None:          
c            CASE (1)
c              Reference plasma:
c            CASE (2)
c             PIN:          
c              source(ic1:ic2) = DBLE(pin(ic1:ic2,ion)%qe) 
c            CASE (3)
c              Calculated...:
            CASEDEFAULT
              CALL User_VolumeEneRecSource(target,source)
          ENDSELECT         
          enerec(ic1:ic2,ion) = source(ic1:ic2)

c...      Ionisation:
          source = 0.0D0
          SELECTCASE (opt%te_ion(target))
            CASE (1000)
c             Preserve (from a previous iteration):
              source(ic1:ic2) = eneion(ic1:ic2,ion)               
            CASE (0)
c             None:          
            CASE (1)
c             Reference plasma:
              source(ic1:ic2) = DBLE(ref_fluid(ic1:ic2,ion)%eneion)   
            CASE (2)
c             PIN:
              source(ic1:ic2) = DBLE(pin(ic1:ic2,ion)%qe) 
            CASE (3)
c             Prescribed:
              CALL SpecifyDistribution(target,-2,0,2,0.1D0,source) 
              source = -1.0D+04 * source 
c              source = -2.25D+07 * source 
            CASEDEFAULT                                            
              CALL User_VolumeEneIonSource(target,source)          
          ENDSELECT         
          eneion(ic1:ic2,ion) = source(ic1:ic2)

c...      Cross-field:
c         ----------------------------------------------------------------

c...      Anomalous:
c         ----------------------------------------------------------------
          source = 0.0D0
          SELECTCASE (opt%te_ano(target))
            CASE (1000)
c             Preserve (from a previous iteration):
              source(ic1:ic2) = eneano(ic1:ic2)               
            CASE (0)
c             None:          
            CASE (1)
c             Reference plasma:
              source(ic1:ic2) = DBLE(ref_fluid(ic1:ic2,ion)%eneano) 
            CASE (2)
c             None, assigned analytically in EvolveTeProfile:          
            CASE (3)
c             None, assigned analytically in EvolveTeProfile:          
            CASEDEFAULT                                            
              CALL ER('AssignEnergySources','Bad TE_ANO',*99)
c              CALL User_VolumeEneAnoSource(target,source)          
          ENDSELECT         
          eneano(ic1:ic2) = source(ic1:ic2)

c...      User defined source:
c         ----------------------------------------------------------------
          eneusr(ic1:ic2) = 0.0D0

        ENDDO
      ENDDO

      CALL IntegrateSources(3)
c...
      IF (log.GT.1) WRITE(logfp,*) 'ASSIGNING ENERGY SOURCES: DONE'

      RETURN
 99   WRITE(0,*) '  TARGET = ',target
      WRITE(0,*) '  TE_ANO = ',opt%te_ano(target)
      STOP
      END

c
c ====================================================================
c
      SUBROUTINE IntegrateSources(mode)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER, INTENT(IN) :: mode

      INTEGER ion
      REAL*8  net,integral(10),ic1,ic2

c...  Additive, so needs to be reset each time:
      SELECTCASE (mode)
c       ------------------------------------------------------------------
        CASE (1)  ! Particles

          DO ion = 1, nion
            IF (iontype(ion).NE.ITY_FLUID) CYCLE

c...        Calculate net particle source along the flux tube:
            parsrc(0      ,ion) = isat(ictarg(LO),ion)
            parsrc(icmax+1,ion) = isat(ictarg(HI),ion)
            parsrc(1:icmax,ion) = parrec(1:icmax,ion) + 
     .                            parion(1:icmax,ion) + 
     .                            parano(1:icmax,ion) 
c...        Integral:
            CALL IntegrateArray(FULL,parsrc(1,ion),1,parint(0,ion))

c...        Detailed check:
            IF (log.GT.0) THEN
              CALL IntegrateArray(FULL,parrec(1,ion),0,integral(1))
              CALL IntegrateArray(FULL,parion(1,ion),0,integral(2))
              CALL IntegrateArray(FULL,parano(1,ion),0,integral(3))
              net = isat(ictarg(LO),ion) + isat(ictarg(HI),ion) + 
     .              integral(1) + integral(2) + integral(3)
              WRITE(logfp,*) 'BIG PARTICLE CHECK:',net,ion
              WRITE(logfp,*) '             isat :',isat(ictarg(LO),ion),
     .                                             isat(ictarg(HI),ion)
              WRITE(logfp,*) '                  :',
     .          (isat(ictarg(LO),ion)+isat(ictarg(HI),ion)),
     .          (isat(ictarg(LO),ion)+isat(ictarg(HI),ion))*ECH
              WRITE(logfp,*) '             rec  :',integral(1)
              WRITE(logfp,*) '             ion  :',integral(2)
              WRITE(logfp,*) '             anp  :',integral(3)

              WRITE(logfp,'(A,1P,E10.2,0P,A)') ' PAR BALANCE:',
     .          (parsrc(0,ion)+parsrc(icmax+1,ion)+parint(TOTAL,ion)) /
     .          (-parsrc(0,ion)-parsrc(icmax+1,ion)) * 100.0D0, ' %'
            ENDIF

          ENDDO
c       ----------------------------------------------------------------
        CASE(2)  ! Momentum

          DO ion = 1, nion
            IF (iontype(ion).NE.ITY_FLUID) CYCLE

c...        Calculate net momentum source along the flux tube:
            momsrc(0      ,ion) = pe(ictarg(LO)) + pi(ictarg(LO),ion)
            momsrc(icmax+1,ion) = pe(ictarg(HI)) + pi(ictarg(HI),ion)
            momsrc(1:icmax,ion) = momvol(1:icmax,ion) + 
     .                            momano(1:icmax,ion)
c...        Integral:
            CALL IntegrateArray(FULL,momsrc(1,ion),1,momint(0,ion))

            IF (log.GT.0) THEN
             WRITE(logfp,*) 'NET MOM CHECK: P0,MAX  :',
     .         momsrc(0      ,ion),
     .         momsrc(icmax+1,ion)
             WRITE(logfp,*) ' NE:',ne(ictarg(LO)    ),ne(ictarg(HI)    )
             WRITE(logfp,*) ' NI:',ni(ictarg(LO),ion),ni(ictarg(HI),ion)
             WRITE(logfp,*) ' VI:',vi(ictarg(LO),ion),vi(ictarg(HI),ion)
             WRITE(logfp,*) ' M :',machno(ictarg(LO),ion),
     .                             machno(ictarg(HI),ion)
             WRITE(logfp,*) ' PE:',pe(ictarg(LO)    ),pe(ictarg(HI)    )
             WRITE(logfp,*) ' PI:',pi(ictarg(LO),ion),pi(ictarg(HI),ion)
             WRITE(logfp,*) ' TE:',te(ictarg(LO)    ),te(ictarg(HI)    )
             WRITE(logfp,*) ' TI:',ti(ictarg(LO),ion),ti(ictarg(HI),ion)
             WRITE(logfp,*) '               INTEGRAL:',momint(TOTAL,ion)
             WRITE(logfp,*) 'RAW VOL:',momvol(1:icmax,ion)
             WRITE(logfp,*) 'RAW ANO:',momano(1:icmax,ion)
             WRITE(logfp,*) 'REF ANO:',ref_fluid(1:icmax,ion)%momano

             WRITE(logfp,'(A,1P,E10.2,0P,A)') '  MOM BALANCE:',
     .        (momsrc(icmax+1,ion) - momsrc(0,ion) - momint(TOTAL,ion))/ 
     .        (momsrc(icmax+1,ion) + momsrc(0,ion)) * 100.0D0,' %'
            ENDIF

          ENDDO
c       ------------------------------------------------------------------
        CASE(3)  ! Energy

          enesrc = 0.0D0
          eneint = 0.0D0

          DO ion = 1, nion
            IF (iontype(ion).NE.ITY_FLUID) CYCLE

c...        Calculate total electron energy source along the flux tube and
c           integrate along the field line:
            enesrc(0      ,1) = enesrc(0      ,1) + 1.0
            enesrc(icmax+1,1) = enesrc(icmax+1,1) + 1.0
            enesrc(1:icmax,1) = enesrc(1:icmax,1) + 
     .                          enerec(1:icmax,ion) +
     .                          eneion(1:icmax,ion) +
     .                          eneano(1:icmax) + 
     .                          eneusr(1:icmax) 

            WRITE(logfp,*) 'ENEION:',eneion(1:icmax,ion)
            WRITE(logfp,*) 'ENEANO:',eneano(1:icmax)
            WRITE(logfp,*) 'ENEUSR:',eneusr(1:icmax)

            CALL IntegrateArray(FULL,enesrc(1,1),1,eneint(0,1))

            ic1 = ictarg(LO)
            ic2 = ictarg(HI)
            qe(ic1) = -eneint(icmid  ,1)
            qe(ic2) = -(eneint(TOTAL,1)-eneint(icmid,1))

            IF (log.GT.0) THEN
              WRITE(logfp,*) 'NET ENE CHECK -TARGETS :',qe(ic1),
     .                                                  qe(ic2)
              WRITE(logfp,*) '              -INTEGRAL:',
     .                eneint(icmid  ,1),
     .                eneint(TOTAL,1)-eneint(icmid,1)
              WRITE(logfp,*) '              -TAR+INT :',
     .                qe(0      )+eneint(icmid  ,1),
     .                qe(icmax+1)+(eneint(TOTAL,1)-eneint(icmid,1))

              WRITE(logfp,*) '              -TOT INT :',eneint(TOTAL,1)
              WRITE(logfp,*) '    SRC:',enesrc(1:icmax,1)
            ENDIF

          ENDDO

        CASE DEFAULT
          CALL ER('IntegrateSources','Unrecognised MODE',*99)
      ENDSELECT

      RETURN
 99   STOP
      END

