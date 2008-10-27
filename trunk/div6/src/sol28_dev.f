c     -*-Fortran-*-
c
c ====================================================================
c
      SUBROUTINE RoutineInDevelopment
      IMPLICIT none
      RETURN
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE LocalGridRefinement
      USE mod_geometry
      IMPLICIT none

      TYPE(type_srf   ) newsrf
      TYPE(type_object) newobj




      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE UpdateTemperatureProfiles
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER ion,ic1,ic2,target,i1,i2,count
      LOGICAL cont
      REAL*8  te1,te2

      cont  = .TRUE.
      count = 0

c...  Reference Te values set in ConserveMomentum (at the moment):
      node(1    )%te = SNGL(te(ictarg(LO)))  ! Need to make more flexible so that it isn't necessarily the target
      node(nnode)%te = SNGL(te(ictarg(HI)))  ! that is floating...

      te1 = te(ictarg(LO))
      te2 = te(ictarg(HI))
 
      i1 = ictarg(LO)
      i2 = ictarg(HI)


      DO WHILE(cont)
        count = count + 1

        WRITE(0,*) 

c...    Te:
        CALL InterpolateTeProfile(2)


c...    Some logic:

        WRITE(0,*) '----------------'
        WRITE(0,*) 'COUNT: ',count
        WRITE(0,*) 'Te:',te(0),te(icmax+1)
        WRITE(0,*) '  :',node(1)%te,node(nnode)%te

        IF (te(i1).LT.te1.AND.te(i2).GT.te2) THEN  ! LEFT OFF -- ADD A DELTA FOR TOLERANCE... 
          WRITE(0,*) 'WHOA!',node(mnode)%icell,node(mnode)%te

          node(mnode)%icell = node(mnode)%icell - 3          
          node(mnode)%te = node(mnode)%te + 0.25

        ELSE

        ENDIF










c...    Ti:   ! Need a separate routine that does this, to be called from various places...
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
        IF (count.GT.2.) THEN
          cont = .FALSE.
        ENDIF

      ENDDO


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE EvolveTeProfile(inode1,inode2,s,target)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none
 
      INTEGER, INTENT(IN) :: inode1,inode2,target
      REAL*8 , INTENT(IN) :: s(0:icmax+1)

      INTEGER :: ic,ic1,ic2,icstep,count,ion
      REAL*8  :: te1,te2,Psol,k,L,tgrad,adjust,vsign,x,frac,tarsign,
     .           flux,qano(0:S28_MAXNKS+1),
     .           qe_src(0:S28_MAXNKS+1),
     .           stotal
      LOGICAL :: cont

      te1 = DBLE(node(inode1)%te)
      te2 = DBLE(node(inode2)%te)
      ic1 = node(inode1)%icell
      ic2 = node(inode2)%icell
      icstep = SIGN(1,ic1-ic2)
      tarsign = DBLE(icstep)

c      IF (opt%bc(target).GE.2.AND.
c     .    (ic1.EQ.ictarg(target).AND.te(ic1).EQ.0.0D0)) THEN
c        WRITE(0,*) 'CYCLING -------------------------'
c        RETURN
c      ENDIF


c      WRITE(0,*) 'TE: ',ic1,ic2,icstep,target
c      WRITE(0,*) '    ',te1,te2

      ion = 1
      k = DBLE(opt%te_kappa(target))
      L = 0.5D0 * DBLE(tube%smax)

c...  Set total power into flux tube:
      SELECTCASE (opt%te_ano_psol(target)) 
        CASE (1000)           ! ie ENEANO array from ref_plasma...? 
          Psol = 0.0D0  ! ???
c          STOP 'OKAY, READY'
        CASE (0)
c...      Conduction model estimate, all in upstream (initial guess only):
          Psol = 2.0D0/7.0D0 * (te2**3.5D0 - te1**3.5D0) * k / L   ! Stangeby, p 190
c          WRITE(0,*) 'PSOL:',psol,ic1,ic2
        CASE (1)
c...      From reference solution:
          Psol = ref_tube%eneano(target)
        CASE DEFAULT
          STOP 'USER OPTION NOT READY, YOU NINNY'
      ENDSELECT

c...  Integrate e/i cf/vol sources/sinks:
      qe_src = 0.0D0
      IF (icstep.EQ.-1) THEN
        qe_src(ic1:ic2) = eneint(ic2,1) - eneint(ic1:ic2,1)  ! Should I set the INT arrays from 0:icmax+1, and have 
        IF (ic1.EQ.0) qe_src(ic1) = eneint(ic2,1)            ! TOTAL = icmax+1 to make things like this simpler?
c        WRITE(0,*) 'QE:',qe(0:icmax+1)
c        WRITE(0,*) 'QE:',enesrc(0:icmax+1,1)
      ELSE
        qe_src(ic2:ic1) = eneint(ic2:ic1,1) - eneint(ic2,1)  
        IF (ic1.EQ.icmax+1) qe_src(ic1) = eneint(TOTAL,1)-eneint(ic2,1)            
c        WRITE(0,*) 'QE:',qe(0:icmax+1)
c        WRITE(0,*) 'QE:',eneint(ic2,1),ic1,icmax+1
c        WRITE(0,*) 'QE:',enesrc(0:icmax+1,1)
      ENDIF

c...  Initializations:
      adjust = 0.0D0
      count = 0
      x = 5.0D0 / 2.0D0
      qano  = 0.0D0
      te(ic2) = te2

c...  Main iteration loop:
      cont = .TRUE.
      DO WHILE (cont)
        cont = .FALSE.
        count = count + 1
 
c...    Setup anomalous energy input:    NEED TO ADD EFFICIENCY, ie DON'T ALWAYS HAVE TO CALCULATE EVERYTHING...
        SELECTCASE (opt%te_ano(target))
          CASE(1000)
c...        Preserved, set in SpecifyEnergySources:
c            qano(ic) = eneano(MIN(ic1,ic2):MAX(ic1,ic2))                      
          CASE(0)
c...        None:
c          CASE(1)
c...        From reference solution:
 !  - does this work here?  need to save QANO into ENEANO... (AND THEN pAss ENEANO BACK...)

          CASE(2)
c...        All power in at the symmetry point:
            qano = Psol 
          CASE(3)
c...        Power distributed uniformly about the symmetry point(!):
            DO ic = ic2+icstep, ic1, icstep  ! No power assigned to the symmetry point cell...
              frac = (s(ic ) - (s(ic2) - 0.5D0 * sdelta(ic2))) / 
     .               (s(ic1) - (s(ic2) - 0.5D0 * sdelta(ic2)))
              qano(ic) = Psol * frac                      
            ENDDO
c          CASE(3)
c...        Between x-points:
c          CASE(4)
c...        B-field depenence... use CalcProfile?   What if heat is not centred at symmetry point? 
          CASE DEFAULT
            STOP 'USER TE OPTION NOT READY'
        ENDSELECT

c...    Electron convected heat flux:
        SELECTCASE (opt%te_conv(target))
          CASE(0)
c...        None:
            DO ic = ic2, ic1, icstep
              qconv(ic) = 0.0D0
            ENDDO
          CASE(1)
            te(ic1) = te1
            DO ic = ic2, ic1, icstep
              IF (te(ic).NE.0.0D0) THEN
                flux = ne(ic) * vi(ic,ion) * tarsign
                qconv(ic) = 5.0D0 / 2.0D0 * ECH * te(ic) * flux
              ENDIF
            ENDDO
          CASE DEFAULT
            STOP 'USER TE OPTION NOT READY'
        ENDSELECT

c...    Ion convected heat flux:


c...    Electron-ion equilibration:
        IF (.FALSE.) THEN
        ENDIF

c...    Calculate Te profile:
        IF (log.GT.0) WRITE(logfp,'(A)') 'TE PROFILE:'
        te(ic1) = 0.0D0
        DO ic = ic2+icstep, ic1, icstep
          qcond(ic) = qano(ic) + qe_src(ic) - qconv(ic) 
          tgrad = -qcond(ic) / (k * te(ic-icstep)**x)                ! Reference...
c...      Note: This is not strictly correct since the convection is a cell centered quantity
c         but the temperature evolution occurs between cell centres, i.e. the convecion
c         term should be a wieghted average, or some such, of the convection
c         across the 2 cells -- similar inconcistencies can likely be found elsewhere
c         and will likely need to be resolved before a full cross-field transport model
c         can work... might even run into trouble before then.  For a similar slop, the
c         treatment of the last half-cell on each ring also needs careful review.
          te(ic) = te(ic-icstep) + (s(ic-icstep) - s(ic)) * tgrad

          qcond(ic) = qcond(ic) * tarsign
          qconv(ic) = qconv(ic) * tarsign

c          qe(ic) = qcond(ic) + qconv(ic)


          IF (log.GE.2) THEN
            WRITE(logfp,'(A,I4,F8.2,1P,2X,2D12.4,2X,2D12.4,1P)')
     .        '  TE_PRO:',ic,
     .        te(ic),
     .        qcond(ic),qconv(ic),
     .        qano(ic),qe_src(ic)
          ENDIF


c          IF (target.EQ.LO)
c     .      WRITE(0,'(A,2I6,2F10.2,1P,6E10.2,2X,2E10.2,0P)') 
c     .        '  Te->',ic,ic-icstep,s(ic),te(ic),vi(ic,ion),flux,
c     .        qano(ic),qe(ic),qconv(ic),qcond(ic),vi(ic,ion),ne(ic)

          IF (te(ic).LE.0.5D0*te1) EXIT  ! There is no way to avoid the near target dip in 
c          IF (te(ic).LT.te1) EXIT       ! Te if Qe is ill-posed?
        ENDDO

c...    Analyse electron profiles and modify parameters, if necessary:
        SELECTCASE (opt%te_ano_psol(target))
          CASE (1000)
          CASE (0)
            IF (DABS(te(ic1)-te1).GT.MIN(0.001D0,0.05D0*te1)) THEN
              vsign = SIGN(1.0D0,te(ic1)-te1)
              IF (adjust.EQ.0.0D0) THEN
                adjust = 0.3D0 * vsign
              ELSEIF ((adjust.GT.0.0D0.AND.vsign.EQ.-1.0D0).OR.
     .                (adjust.LT.0.0D0.AND.vsign.EQ. 1.0D0)) THEN
                adjust = -0.3D0 * adjust
              ENDIF
c...          Scale Psol to improve match to specified target temperature:
              SELECTCASE (1)
                CASE (1)  ! Power into SOL:
                  Psol = Psol * (1.0D0 + adjust)
                CASE DEFAULT
                  STOP 'NOT READY, SORRY'
              ENDSELECT
              cont = .TRUE.
              IF (log.GE.2) THEN
                WRITE(logfp,'(A,3F11.3,1P,2E10.2,0P)') 
     .            'ADJUST:',te(ic1),te1,vsign,adjust,Psol
              ENDIF
            ENDIF
          CASE DEFAULT
            CALL ER('EvolveTeProfile','Bad TE_ANO_PSOL',*99)
        ENDSELECT
c...    Analyse ion profiles:



c...    If convection being used, need to run for extra iterations
c       to make sure the solution is converged:


      ENDDO ! End of main iterative loop


      WRITE(0,*) 'QCOND:',qcond(ic1),qcond(ic2)
      WRITE(0,*) 'TE   :',te(ic1),te(ic2)
      WRITE(0,*) 'PSOL :',psol,ic1,ic2,opt%te_ano_psol(target)

c...  Store Psol information:
      tube%eneano(target) = Psol   ! Keep this up...?

c...  Linear distribution of anomalous power for now...
      SELECTCASE (opt%te_ano(target))
        CASE(1000)
        CASE(2)
          eneano(ic2+icstep) = Psol / sdelta(ic2+icstep)  ! Approximate..!
          WRITE(0,*) 'ENEA:',psol,ic2+icstep
        CASE(3)
          IF (ic1.GT.ic2) THEN
            ic  = ic1
            ic1 = ic2
            ic2 = ic
          ENDIF
          ic1 = MAX(ic1,1)
          ic2 = MIN(ic2,icmax)
          IF (target.EQ.LO) ic2 = ic2 - 1
          IF (target.EQ.HI) ic1 = ic1 + 1
          stotal = SUM(sdelta(ic1:ic2))
          eneano(ic1:ic2) = Psol / stotal
          WRITE(0,*) 'ENEA:',psol,ic1,ic2,icstep,stotal
        CASE DEFAULT
          CALL ER('EvolveTeProfile','Bad TE_ANO',*99)
      ENDSELECT



      RETURN
 99   STOP
      END





