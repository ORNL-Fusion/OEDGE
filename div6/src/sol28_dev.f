c     -*Former Mode Specification*-
c
c ====================================================================
c
      SUBROUTINE LineCutTube(p1,p2,itube,p)
      USE mod_geometry
      USE mod_sol28_global
      IMPLICIT none

      INTEGER, INTENT(IN)  :: itube      
      REAL*8 , INTENT(IN)  :: p1(2),p2(2)
      REAL*8 , INTENT(OUT) :: p(2)

      INTEGER GetObject

      INTEGER ic1,ic2,icell,iobj
      REAL*8  p3(2),p4(2),t12,t34

      ic1 = tube(itube)%cell_index(1)
      ic2 = tube(itube)%cell_index(2)
      DO icell = ic1, ic2
        iobj = GetObject(icell,IND_CELL)                   
        CALL GetVertex(iobj,2,p3(1),p3(2))          
        CALL GetVertex(iobj,3,p4(1),p4(2))          
        CALL CalcInter(p1(1),p1(2),p2(1),p2(2),
     .                 p3(1),p3(2),p4(1),p4(2),t12,t34)
        IF (t12.GE.0.0D0.AND.t34.GE.0.0D0.AND.t34.LT.1.0D0) THEN
          p(1) = p3(1) + t34 * (p4(1) - p3(1))
          p(2) = p3(2) + t34 * (p4(2) - p3(2))
          EXIT
        ENDIF
      ENDDO
      IF (icell.EQ.ic2+1) 
     .  CALL ER('LineCutTube','Intersection not found',*99)

      RETURN
 99   WRITE(0,*) '  ITUBE = ',itube
      STOP
      END
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
      SUBROUTINE CalculateTemperatureProfile(inode1,inode2,s,target,
     .                                       mode,t)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none
 
      INTEGER, INTENT(IN ) :: inode1,inode2,target,mode
      REAL*8 , INTENT(IN ) :: s(0:icmax+1)
      REAL*8 , INTENT(OUT) :: t(*)

      INTEGER ion,ic,ic1,ic2,i1,i2,count,adjust_count
      INTEGER icstep
      LOGICAL cont
      REAL*8  te1,te2,Psol,L,k,stotal,vsign,adjust,dist,Pano

      ion = 1

      ic1 = node(inode1)%icell
      ic2 = node(inode2)%icell
      icstep = SIGN(1,ic1-ic2)
c      te1 = DBLE(node(inode1)%te)

      IF (mode.EQ.1) THEN 
        IF (opt%bc(target).EQ.3) THEN
          IF (target.EQ.LO) te1 = te(ictarg(LO))
          IF (target.EQ.HI) te1 = te(ictarg(HI))
        ELSE
          te1 = DBLE(node(inode1)%te)
        ENDIF
        te2 = DBLE(node(inode2)%te)
c        IF (opt%bc(target).EQ.3) THEN
c          IF (target.EQ.LO) te2 = te(ictarg(LO))
c          IF (target.EQ.HI) te2 = te(ictarg(HI))
c        ELSE
c          te2 = DBLE(node(inode2)%te)
c        ENDIF
        k = DBLE(opt%te_kappa(target))
      ELSEIF (mode.EQ.2) THEN 
        IF (opt%bc(target).EQ.3) THEN
          IF (target.EQ.LO) te1 = ti(ictarg(LO),ion)
          IF (target.EQ.HI) te1 = ti(ictarg(HI),ion)
        ELSE
          te1 = DBLE(node(inode1)%ti(ion))
        ENDIF
        te2 = DBLE(node(inode2)%ti(ion))
        k = DBLE(opt%ti_kappa(target))
      ELSE
        STOP 'CRAP MAN COMBO!'
      ENDIF

      L = 0.5D0 * DBLE(tube%smax)

      adjust_count = 0
      adjust = 0.0D0
      count  = 0
      dist   = 1.0D0

      WRITE(logfp,'(A,I6)') 'TARGET:',target

      cont = .TRUE.
      DO WHILE (cont)
        cont = .FALSE.

        count = count + 1

c...LEFT OFF -- NEED TO PUT THIS IN A SEPARATR SUBROUTINE, AND TO GO THORUHG AND CLEAN UP THE BLOODY CODE
c   STO THAT I KNWO THINGS ARE DORN RIGHT!

c...    Set total power into flux tube:
        SELECTCASE (opt%te_ano_psol(target)) 
          CASE (1000)           ! ie ENEANO array from ref_plasma...? 
            Psol = 0.0D0  ! ???
            Pano = Psol
c            STOP 'OKAY, READY'
          CASE (0)
c...        Conduction model estimate, all in upstream (initial guess only):
            Psol = 2.0D0/7.0D0 * (te2**3.5D0 - te1**3.5D0) * k / L   ! Stangeby, p 190
            Pano = Psol
c            WRITE(0,*) 'PSOL:',psol,ic1,ic2
          CASE (1)
c...        From reference solution:
            Psol = ref_tube%eneano(target)
            Pano = Psol
          CASE (2)
c...        From conservation:
            Psol = 8.0D0 * te1 * 
     .             DABS(isat(ictarg(target),ion)) * ECH
            write(logfp,*) 'target heat flux',target,psol
            IF (target.EQ.LO) THEN 
              Psol = Psol - dist * eneint(icmid  ,1)
            ELSE
              Psol = Psol - dist * (eneint(TOTAL,1) - eneint(icmid,1))
            ENDIF
            Pano = Psol
          CASE (3)
c...        From conservation:
            IF (target.EQ.LO) THEN 
              Psol = -eneint(icmid,1)
            ELSE
              Psol = -(eneint(TOTAL,1) - eneint(icmid,1))
            ENDIF

            Pano = 8.0D0 * te1 * 
     .              DABS(isat(ictarg(target),ion)) * ECH

            write(logfp,*) 'target heat flux',target,pano

            Pano = Pano + Psol * (1.0D0 + dist)
          CASE DEFAULT
            STOP 'USER OPTION NOT READY, YOU NINNY'
        ENDSELECT

        CALL EvolveTemperatureProfile(inode1,inode2,s,target,Psol,Pano,
     .                                te1,te2,dist,mode,t)

c        if (count.EQ.30) stop 'check 1'

c...    Analyse electron profiles and modify parameters, if necessary:

        IF (DABS(t(ic1)-te1).GT.MIN(1.0D0,0.05D0*te1)) THEN
c        IF (DABS(t(ic1)-te1).GT.MIN(0.001D0,0.05D0*te1)) THEN
c          vsign = SIGN(1.0D0,t(ic1)-te1) 
c          IF (target.EQ.LO) vsign = vsign * SIGN(1.0D0,te2-te1)
          vsign = SIGN(1.0D0,t(ic1)-te1) * SIGN(1.0D0,te2-te1)
          IF     (adjust.EQ.0.0D0) THEN
            adjust = 0.3D0 * vsign  
          ELSEIF ((adjust.GT.0.0D0.AND.vsign.EQ.-1.0D0).OR.   
     .            (adjust.LT.0.0D0.AND.vsign.EQ. 1.0D0)) THEN
            adjust = -0.3D0 * adjust
            adjust_count = -1
          ELSE
            IF (adjust_count.GE.0 ) adjust_count = adjust_count + 1
            IF (adjust_count.EQ.10) THEN  ! Only "kick" if ADJUST is not oscillating
              adjust = adjust * 10.0D0
              adjust_count = 0
              WRITE(logfp,*) '*** adjust kick! ***'
            ENDIF
          ENDIF
          cont = .TRUE.
        ENDIF

c...    Analyse electron profiles and modify parameters, if necessary:

        SELECTCASE (opt%bc(target))  ! (opt%te_ano_psol(target))  ! *** CHANGE VARIABLE *** 
          CASE (1)
c...        Adjust the anomalous power term on the half ring:
c...        Scale Pano to improve match to specified target temperature:
            SELECTCASE (1)
              CASE (1)  ! Power into SOL:
c                dist = dist * (1.0D0 - adjust)
                dist = dist + adjust  ! new
c                dist = dist * (1.0D0 + adjust)  ! old
c                IF (dist.LT.0.0D0.OR.dist.GT.1.0D0) THEN
c                  WRITE(0,*) 'KARLS NEW SOLVER IS CRAP'
c                  STOP
c                ENDIF
c                Pano = Pano * (1.0D0 + adjust)
              CASE DEFAULT
                STOP 'NOT READY, SORRY'
            ENDSELECT
            IF (logop.GE.2) THEN
              WRITE(logfp,'(A,I3,5F11.4,1P,2E14.6,0P,F12.6)') 
     .          'ADJUST:',count,t(ic1),te1,DABS(t(ic1)-te1),
     .                 MIN(1.0D0,0.05D0*te1),vsign,adjust,Pano,dist
            ENDIF

          CASE (3)

            node(inode2)%te = node(inode2)%te - adjust
            te2 = DBLE(node(inode2)%te)

            IF (logop.GE.2) THEN
              WRITE(logfp,'(A,3F11.4,1P,E14.6,0P,F11.4)') 
     .          'ADJUST:',t(ic1),te1,vsign,adjust,te2
            ENDIF

          CASE DEFAULT
            CALL ER('EvolveTeProfile','Bad BC',*99)
        ENDSELECT

c...    Analyse ion profiles:

c...    If convection being used, need to run for extra iterations
c       to make sure the solution is converged:

        IF (count.EQ.150) THEN
          t(ic1:ic2) = te1
          WRITE(0,*) '*** ENERGY MODEL FAILED ***'
          WRITE(logfp,*) '*** ENERGY MODEL FAILED ***'
          EXIT
c          STOP 'EXCESS ITERATIONS'
        ENDIF

      ENDDO


c...  Store Pano information:
      tube%eneano(target) = Pano   ! Keep this up...?

c...  Linear distribution of anomalous power for now...
      SELECTCASE (opt%te_ano(target))
        CASE(1000)
        CASE(1)
        CASE(2)
c...      Everything in at the midplane:
          eneano(ic2+icstep) = Pano / sdelta(ic2+icstep)  ! Approximate..!
c          WRITE(0,*) 'ENEA:',pano,ic2+icstep
        CASE(3,5)
c...      Distributed uniformly along the half ring:
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
          eneano(ic1:ic2) = Pano / stotal
c          WRITE(0,*) 'ENEA:',pano,ic1,ic2,icstep,stotal
        CASE DEFAULT
          CALL ER('EvolveTeProfile','Bad TE_ANO',*99)
      ENDSELECT


      IF (logop.GE.2) THEN
        WRITE(logfp,'(A)') 'DONE TE ITERATIONS'
      ENDIF

      IF (target.EQ.HI.AND.opt%bc(target).EQ.1) THEN
        WRITE(logfp,*) 'INTEGRATING ENERGY SOURCES AGAIN:'
        CALL IntegrateSources(3)

        ic1 = ictarg(LO)
        ic2 = ictarg(HI)

c        qe(ic1) = -eneint(icmid  ,1)
c        qe(ic2) = -(eneint(TOTAL,1)-eneint(icmid,1))

        WRITE(logfp,*) 'QE   :',qe(0),qe(icmax+1)
        WRITE(logfp,*) 'ISAT :',isat(ictarg(LO:HI),ion)
        WRITE(logfp,*) 'TE   :',te(ictarg(LO:HI))
        WRITE(logfp,*) 'TI   :',ti(ictarg(LO:HI),ion)

        STOP 'need to check all this for Ti, not just Te'
        gamma(LO) = qe(0)       / (isat(ic1,ion)*ECH) / te(ic1)
        gamma(HI) = qe(icmax+1) / (isat(ic2,ion)*ECH) / te(ic2)

c        gamma = 3.0D0
c        qe(ic1) = 3.0D0 * (ECH*te(ic1)) * isat(ic1,ion)
c        qe(ic2) = 3.0D0 * (ECH*te(ic2)) * isat(ic1,ion)

        WRITE(logfp,*) 'GAMMA:',gamma(LO:HI)
        WRITE(logfp,*) '...and again...'

        CALL IntegrateSources(3)
c        gamma(LO:HI) = 3.0D0
c        STOP 'sdgsdg'
      ENDIF


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE EvolveTemperatureProfile(inode1,inode2,s,target,Psol,
     .                                    Pano,te1,te2,dist,mode,t)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none
 
      INTEGER, INTENT(IN ) :: inode1,inode2,target,mode
      REAL*8 , INTENT(IN ) :: s(0:icmax+1),Psol,Pano,te1,te2,dist
      REAL*8 , INTENT(OUT) :: t(*)

      INTEGER :: ic,ic1,ic2,ic3,icstep,count,ion
      REAL*8  :: tgrad,x,k,frac,tarsign,peak_pos,
     .           flux,qano(0:S28_MAXNKS+1),
     .           qe_src(0:S28_MAXNKS+1),
     .           qconv_e(0:S28_MAXNKS+1),qconv_i(0:S28_MAXNKS+1)
      LOGICAL :: cont

c      te1 = DBLE(node(inode1)%te)
c      te2 = DBLE(node(inode2)%te)
      ic1 = node(inode1)%icell
      ic2 = node(inode2)%icell
      icstep = SIGN(1,ic1-ic2)
      tarsign = DBLE(icstep)

c      IF (opt%bc(target).GE.2.AND.
c     .    (ic1.EQ.ictarg(target).AND.te(ic1).EQ.0.0D0)) THEN
c        WRITE(0,*) 'CYCLING -------------------------'
c        RETURN
c      ENDIF


c      WRITE(logfp,*) 'TE: ',ic1,ic2,icstep,target
c      WRITE(logfp,*) '    ',te1,te2

      ion = 1

c...  Integrate e/i cf/vol sources/sinks:
      qe_src = 0.0D0

c...  Set total power into flux tube:
      SELECTCASE (opt%te_ano_psol(target)) 
        CASE (3)
          IF (icstep.EQ.-1) THEN
            qe_src(ic1:ic2) = (eneint(ic2,1) - eneint(ic1:ic2,1))  
            IF (ic1.EQ.0) qe_src(ic1) = eneint(ic2,1)              

            peak_pos = 0.03D0
            DO ic = 1, ic2
              IF ( (cell(ic)%s / tube%smax).GT.peak_pos ) EXIT
            ENDDO
            IF (ic.EQ.ic2+1) 
     .        CALL ER('EvolveTeProfile','Out of range 1',*99)
            
            qe_src(ic1:ic) = qe_src(ic1:ic) - Psol * dist

            write(logfp,*) 'ic,psol',ic,ic1,ic2,psol

c            WRITE(0,*) 'QE:',qe(0:icmax+1)
c            WRITE(0,*) 'QE:',enesrc(0:icmax+1,1)
          ELSE
            qe_src(ic2:ic1) = (eneint(ic2:ic1,1) - eneint(ic2,1))
            IF (ic1.EQ.icmax+1) 
     .        qe_src(ic1) = (eneint(TOTAL,1)-eneint(ic2,1))

            peak_pos = 0.01D0
            DO ic = ic1-1, ic2, -1
c             write(logfp,*) 'frac',cell(ic)%s / tube%smax,1.0D0-peak_pos
              IF ( (cell(ic)%s / tube%smax).LT.(1.0D0-peak_pos) ) EXIT
            ENDDO
            IF (ic.EQ.ic2-1) 
     .        CALL ER('EvolveTeProfile','Out of range 2',*99)
            
            qe_src(ic:ic1) = qe_src(ic:ic1) - Psol * dist

            write(logfp,*) 'ic,psol',ic,ic1,ic2,psol

c            WRITE(0,*) 'QE:',qe(0:icmax+1)
c            WRITE(0,*) 'QE:',eneint(ic2,1),ic1,icmax+1
c            WRITE(0,*) 'QE:',enesrc(0:icmax+1,1)
           ENDIF
        CASE DEFAULT
          IF (icstep.EQ.-1) THEN
            qe_src(ic1:ic2) = dist * (eneint(ic2,1) - eneint(ic1:ic2,1))  ! Should I set the INT arrays from 0:icmax+1, and have 
            IF (ic1.EQ.0) qe_src(ic1) = dist * eneint(ic2,1)              ! TOTAL = icmax+1 to make things like this simpler?
c            WRITE(0,*) 'QE:',qe(0:icmax+1)
c            WRITE(0,*) 'QE:',enesrc(0:icmax+1,1)
          ELSE
            qe_src(ic2:ic1) = dist * (eneint(ic2:ic1,1) - eneint(ic2,1))
            IF (ic1.EQ.icmax+1) 
     .        qe_src(ic1) = dist * (eneint(TOTAL,1)-eneint(ic2,1))
c            WRITE(0,*) 'QE:',qe(0:icmax+1)
c            WRITE(0,*) 'QE:',eneint(ic2,1),ic1,icmax+1
c            WRITE(0,*) 'QE:',enesrc(0:icmax+1,1)
          ENDIF
      ENDSELECT

c...  Initializations:
      IF (mode.EQ.1) THEN
        k = DBLE(opt%te_kappa(target))
      ELSE
        k = DBLE(opt%ti_kappa(target))
      ENDIF
      x = 5.0D0 / 2.0D0
      qano  = 0.0D0
      t(ic2) = te2
 
c...  Setup anomalous energy input:    NEED TO ADD EFFICIENCY, ie DON'T ALWAYS HAVE TO CALCULATE EVERYTHING...
      SELECTCASE (opt%te_ano(target))
        CASE(1000)
        CASE(0)
c...      None:
        CASE(1)
c...      From reference solution - ANO source incorporated into QE_SRC via ENEINT:
        CASE(2)
c...      All power in at the symmetry point:
          qano = Pano 
        CASE(3)
c...      Power distributed uniformly about the symmetry point(!):
          DO ic = ic2+icstep, ic1, icstep  ! No power assigned to the symmetry point cell...
            frac = (s(ic ) - (s(ic2) - 0.5D0 * sdelta(ic2))) / 
     .             (s(ic1) - (s(ic2) - 0.5D0 * sdelta(ic2)))
            qano(ic) = Pano * frac                      
          ENDDO
c        CASE(3)
c...      Between x-points:
c        CASE(4)
c...      B-field depenence... use CalcProfile?   What if heat is not centred at symmetry point? 
        CASE(5)
c...      Power distribution a function of iteration parameter:
          qano = Pano 
          IF (target.EQ.LO) THEN
            ic3 = ic2 - NINT(dist * DBLE(REAL(ic2 - ic1)) )
          ELSE
            ic3 = ic2 + NINT(dist * DBLE(REAL(ic1 - ic2)) )
          ENDIF
          WRITE(logfp,*) 'IC1-3',ic1,ic2,ic3
          DO ic = ic2+icstep, ic3, icstep  ! No power assigned to the symmetry point cell...
            frac = (s(ic ) - (s(ic2) - 0.5D0 * sdelta(ic2))) / 
     .             (s(ic3) - (s(ic2) - 0.5D0 * sdelta(ic2)))
            qano(ic) = Pano * frac                      
          ENDDO
        CASE DEFAULT
          STOP 'USER TE OPTION NOT READY'
      ENDSELECT

c...  Electron convected heat flux:
      IF (mode.EQ.1) THEN 
        SELECTCASE (opt%te_conv(target))
          CASE(0)
c...        None:
            DO ic = ic2, ic1, icstep
              qconv_e(ic) = 0.0D0
            ENDDO
          CASE(1)
            te(ic1) = te1
            DO ic = ic2, ic1, icstep
c              WRITE(logfp,*) 'CONV:',ic,te(ic),ne(ic),vi(ic,ion)
              IF (te(ic).NE.0.0D0) THEN
                flux = ne(ic) * vi(ic,ion) * tarsign
                qconv_e(ic) = 5.0D0 / 2.0D0 * ECH * te(ic) * flux
              ENDIF
            ENDDO
          CASE DEFAULT
            STOP 'USER TE OPTION NOT READY'
        ENDSELECT
      ELSE
        qconv_e = 0.0D0
      ENDIF

c...  Ion convected heat flux:
      IF (mode.EQ.2) THEN 
        SELECTCASE (opt%ti_conv(target))  ! Don't need this when solving the electron channel...
          CASE(0)
c...        None:
            DO ic = ic2, ic1, icstep
              qconv_i(ic) = 0.0D0
            ENDDO
          CASE(1)
            te(ic1) = te1                        ! *** ASSUMES Ti = Te, on this line and 4 below ***
            DO ic = ic2, ic1, icstep
              IF (te(ic).NE.0.0D0) THEN
                flux = ne(ic) * vi(ic,ion) * tarsign
                qconv_i(ic) = (5.0D0 / 2.0D0 * ECH * ti(ic,ion) + 
     .                         0.5D0 * mi(ion) * vi(ic,ion)**2) * flux
              ENDIF
            ENDDO
          CASE DEFAULT
            STOP 'USER TI OPTION NOT READY'
        ENDSELECT
      ELSE
        qconv_i = 0.0D0
      ENDIF

      qconv = qconv_e + qconv_i

c...  Electron-ion equilibration:
      IF (.FALSE.) THEN
      ENDIF

c...  Calculate Te profile:
      t(ic1) = 0.0D0
      DO ic = ic2+icstep, ic1, icstep
        qcond(ic) = qano(ic) + qe_src(ic) - qconv(ic) 
        tgrad = -qcond(ic) / (k * t(ic-icstep)**x)                ! Reference...

c...    Note: This is not strictly correct since the convection is a cell centered quantity
c       but the temperature evolution occurs between cell centres, i.e. the convecion
c       term should be a wieghted average, or some such, of the convection
c       across the 2 cells -- similar inconcistencies can likely be found elsewhere
c       and will likely need to be resolved before a full cross-field transport model
c       can work... might even run into trouble before then.  For a similar slop, the
c       treatment of the last half-cell on each ring also needs careful review.

        t(ic) = t(ic-icstep) + (s(ic-icstep) - s(ic)) * tgrad

        qcond(ic) = qcond(ic) !* tarsign
        qconv(ic) = qconv(ic) !* tarsign

c        qe(ic) = qcond(ic) + qconv(ic)

        IF (logop.GE.2) THEN
          IF (ic.EQ.ic2+icstep) THEN
            WRITE(logfp,'(A4,A8,2X,3A12,2X,2A12,2X,2A12)')
     .        'IND','Te','qcond','qconv_e','qconv_i','qano','qe_src',
     .        'net','pow'
          ENDIF
          WRITE(logfp,'(I4,F8.2,2X,1P,3D12.4,2X,2D12.4,2X,2D12.4,2X,1P,
     .                  2E10.2)')
     .      ic,
     .      t(ic),
     .      qcond(ic),qconv_e(ic),qconv_i(ic),
     .      qano(ic),qe_src(ic),
     .      qcond(ic)+qconv_e(ic)+qconv_i(ic) - 
     .      (qano(ic) + qe_src(ic)),
     .      qano(ic) + qe_src(ic),
     .      tgrad,(s(ic-icstep) - s(ic))
        ENDIF
c          IF (target.EQ.LO)
c     .      WRITE(0,'(A,2I6,2F10.2,1P,6E10.2,2X,2E10.2,0P)') 
c     .        '  Te->',ic,ic-icstep,s(ic),t(ic),vi(ic,ion),flux,
c     .        qano(ic),qe(ic),qconv(ic),qcond(ic),vi(ic,ion),ne(ic)
c        IF (t(ic).LE.0.1D0*te1) EXIT  
c        IF (t(ic).LE.0.5D0*te1.OR.t(ic).GT.1.5D0*te2) THEN
c          t(ic) = MIN(MAX(t(ic),te1),te2)
c          EXIT  ! There is no way to avoid the near target dip in 
c        ENDIF
        IF (t(ic).LE.0.5D0*MIN(te1,te2)) EXIT  ! There is no way to avoid the near target dip in 
c        IF (t(ic).LE.0.5D0*te1) EXIT  ! There is no way to avoid the near target dip in 
c        IF (t(ic).LT.te1) EXIT        ! Te if Qe is ill-posed?
      ENDDO


c      WRITE(0,*) 'QCOND:',qcond(ic1),qcond(ic2)
c      WRITE(0,*) 'TE   :',t(ic1),t(ic2)
c      WRITE(0,*) 'PANO :',pano,ic1,ic2,opt%te_ano_psol(target)



      RETURN
 99   STOP
      END





