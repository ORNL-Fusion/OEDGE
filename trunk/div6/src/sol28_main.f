c     -*Fortran*-
c
c ====================================================================
c
      SUBROUTINE FloatingTargets
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER ion,target,inode,ict1,ict2,ict,off_target
      REAL*8  net,integral(10),p1,p2,source(icmax),pmatch,
     .        srcint(0:icmax),delta_p,isat_total,
     .        alpha(2),beta(2),delta,A,B,C,D,a1,b1,c1,radical,cs
 
      IF (logop.GT.1) WRITE(logfp,*) 'FLOATING TARGETS: BEGIN'

      DO ion = 1, nion
        IF (iontype(ion).NE.ITY_FLUID) CYCLE

c...    Assign target particle flux is required (can't do this before
c       now since the upstream sources have to be assigned for DELTA_P):
c       --------------------------------------------------------------

        DO target = LO, HI
          ict = ictarg(target)

          SELECTCASE (opt%bc(target))
            CASE (1)
c             Do nothing:
            CASE (2)
              CALL IntegrateArray(FULL,parrec(1,ion),0,integral(1))
              CALL IntegrateArray(FULL,parion(1,ion),0,integral(2))
              CALL IntegrateArray(FULL,parano(1,ion),0,integral(3))

              ict1 = ictarg(LO)
              ict2 = ictarg(HI)

              IF (target.EQ.LO) off_target = HI
              IF (target.EQ.HI) off_target = LO

              IF (target.EQ.LO.AND.opt%bc(HI).EQ.1.OR.
     .            target.EQ.HI.AND.opt%bc(LO).EQ.1) THEN

                isat_total = integral(1) + integral(2) + integral(3) + 
     .                       isat(off_target,ion) 

                isat(ict,ion) = -isat_total  ! Okay...?
              ELSE
c...            Both targets are floating - just jsat, soon to be obsolete:

                CALL IntegrateArray(FULL,momvol(1,ion),0,integral(4))
                CALL IntegrateArray(FULL,momano(1,ion),0,integral(5))

                IF (target.EQ.LO) delta_p = -integral(4) - integral(5)
                IF (target.EQ.HI) delta_p =  integral(4) + integral(5)

                alpha(target) = delta_p * ECH / DSQRT(mi(ion))

                beta(LO) = DSQRT(ECH * (te(ict1) + ti(ict1,ion))) * 
     .                     (1.0D0 + machno(ict1,ion)**2) /
     .                     machno(ict1,ion)
                beta(HI) = DSQRT(ECH * (te(ict2) + ti(ict2,ion))) *
     .                     (1.0D0 + machno(ict2,ion)**2) /
     .                     machno(ict2,ion)

                isat_total = (integral(1) + integral(2) + 
     .                        integral(3)) * ECH

                isat(ict,ion) = -1.0D0 *  
     .            (alpha(target) + 
     .             isat_total * beta(off_target)) /
     .            (beta(LO) + beta(HI)) / ECH
              ENDIF

              ne(ict)     = DABS(isat(ict,ion) / vi(ict,ion))    ! These assignments should match the 
              ni(ict,ion) = ne(ict)                              ! ones in AssignTargetValues in main.f --
              pe(ict)     = ne(ict) * te(ict) * ECH              ! they should really be put into a 
              pi(ict,ion) = ni(ict,ion) *                        ! separate routine and called from here and there...
     .                     (ti(ict,ion)*ECH + mi(ion) * vi(ict,ion)**2)

            CASE (3)
              CALL IntegrateArray(FULL,parsrc(1,ion),0,integral(1))
c              CALL IntegrateArray(FULL,parion(1,ion),0,integral(2))
c              CALL IntegrateArray(FULL,parano(1,ion),0,integral(3))

              ict1 = ictarg(LO)
              ict2 = ictarg(HI)

              IF (target.EQ.LO) off_target = HI
              IF (target.EQ.HI) off_target = LO

              IF (target.EQ.LO.AND.opt%bc(HI).EQ.1.OR.
     .            target.EQ.HI.AND.opt%bc(LO).EQ.1) THEN
                STOP 'THIS CODE NOT READY'
              ELSE
c...            Both targets are floating:

                CALL IntegrateArray(FULL,momsrc(1,ion),0,integral(4))
c                CALL IntegrateArray(FULL,momvol(1,ion),0,integral(4))
c                CALL IntegrateArray(FULL,momano(1,ion),0,integral(5))

                IF (target.EQ.LO) delta_p = -integral(4)
                IF (target.EQ.HI) delta_p =  integral(4)
c                IF (target.EQ.LO) delta_p = -integral(4) - integral(5)
c                IF (target.EQ.HI) delta_p =  integral(4) + integral(5)

                delta = delta_p * DSQRT(ECH / mi(ion))

                alpha(LO:HI) = 2.0D0  ! = (1.0 + Ti/Te)

                beta(LO) = DSQRT(DABS(qe(ict1)) * alpha(LO)/gamma(LO)) *
     .                     (1.0D0 + machno(ict1,ion)**2) /
     .                     machno(ict1,ion)
                beta(HI) = DSQRT(DABS(qe(ict2)) * alpha(HI)/gamma(HI)) *
     .                     (1.0D0 + machno(ict2,ion)**2) /
     .                     machno(ict2,ion)

                isat_total = integral(1) * ECH
c                isat_total = (integral(1) + integral(2) + 
c    .                        integral(3)) * ECH

                A = beta(target)**2 - beta(off_target)**2
                B = -2.0D0 * beta(target) * beta(off_target) 
                C = isat_total
                D = beta(off_target)**2 * C - delta**2

                a1 = A**2 + B**2
                b1 = 2.0D0 * A * D - B**2 * C
                c1 = D**2

                radical = b1**2 - 4.0D0 * a1 * c1


                WRITE(logfp,*) 
                WRITE(logfp,'(1X,A,I2)') 'DATA ',target
                WRITE(logfp,*) 'gamma   :',gamma(LO:HI)          
                WRITE(logfp,*) 'ict     :',ict,target,ion          
                WRITE(logfp,*) 'parsrc  :',integral(1)  
                WRITE(logfp,*) 'isat_t  :',isat_total              
                WRITE(logfp,'(A,1P,3E12.4,0P)') ' P      :',delta_p,
     .                                      integral(4:5)
                WRITE(logfp,*) 'beta    :',beta(target),beta(off_target)
                WRITE(logfp,*) 'delta   :',delta
                WRITE(logfp,*) 'machno  :',machno(ict1,ion),
     .                                     machno(ict2,ion)
                WRITE(logfp,*) 'te      :',te(ict)
                WRITE(logfp,*) 'qe      :',qe(ict1),qe(ict2)
                WRITE(logfp,'(A,1P,4E12.4,0P)') 'A,B,C,D :',A,B,C,D
                WRITE(logfp,'(A,1P,3E12.4,0P)') 'A1,B1,C1:',a1,b1,c1
                WRITE(logfp,*) 'RADICAL :',radical
                WRITE(logfp,*) 'D/A     :',D/A

                IF (radical.GE.0.0D0) THEN

                  isat(ict,ion) = (b1 - DSIGN(1.0D0,delta_p) * 
     .                             DSQRT(radical)) / (2.0D0 * a1 * ECH)

                  WRITE(logfp,*) 'isat    :',isat(ict,ion) * ECH
                  WRITE(logfp,*) 'te      :',
     .              DABS(qe(ict)/isat(ict,ion))/ECH/gamma(target)
                  WRITE(logfp,*)

                ELSE
                  STOP 'TROUBLES MAN'
                ENDIF

                IF (target.EQ.HI) THEN 
                  WRITE(logfp,*) 
     .              'isat tot:',SUM(isat(ictarg(LO:HI),ion))
                  WRITE(logfp,*) 
                ENDIF

              ENDIF

              te(ict)     = DABS(qe(ict) / (isat(ict,ion) * ECH)) / 
     .                      gamma(target)
              ti(ict,ion) = (alpha(target) - 1.0D0) * te(ict)

              cs = DSQRT( (te(ict) + ti(ict,ion)) * ECH / mi(ion) )  ! Needs improvement... CalcCs
              vi(ict,ion) = cs * machno(ict,ion) * tsign(target)

              ne(ict)     = DABS(isat(ict,ion) / vi(ict,ion))      ! These assignments should match the 
              ni(ict,ion) = ne(ict)                                ! ones in AssignTargetValues in main.f --
              pe(ict)     = ne(ict) * te(ict) * ECH                ! they should really be put into a 
              pi(ict,ion) = ni(ict,ion) *                          ! separate routine and called from here and there...
     .                     (ti(ict,ion)*ECH + mi(ion) * vi(ict,ion)**2)
             
            CASE DEFAULT
              CALL ER('ConserveParticles','Bad target option',*99)
          ENDSELECT

        ENDDO

      ENDDO

      IF (logop.GT.1) WRITE(logfp,*) 'FLOATING TARGETS: DONE'
      
      RETURN
 99   STOP
      END
c
c ====================================================================
c
c subroutine: InterpolateReferencePlasma
c
c  Treats the core, SOL and PFR separately, i.e. the interpolation 
c  will make a step at the separatrix if the focus grid is narrower 
c  there than the reference grid/plasma.
c 
      SUBROUTINE InterpolateReferencePlasma(tube,nion,fluid,cell)
     .                                      
      USE mod_sol28_params
      USE mod_sol28_reference
c      USE mod_sol28_global
      USE mod_solps_params
      USE mod_solps
      IMPLICIT none

      INTEGER          :: nion
      TYPE(type_tube ) :: tube
      TYPE(type_fluid) :: fluid(tube%n,nion)
      TYPE(type_cell ) :: cell (tube%n)

      REAL, PARAMETER :: TOL = 1.0E-03 , ECH = 1.6022E-19, 
     .                   AMU = 1.67E-27

      INTEGER ion,it,it1,it2,ref_ic,ref_ic1,ref_ic2,int_ic1,int_ic2,
     .        i1,ic,ic1,ic2,itarget,n,i
      LOGICAL output
      REAL    pfr,fr,
     .        val1(5,1000),val2(5,1000),val3(5,1000),val4(5,1000),
     .        val5(5,2),val6(5,2),
     .        mi,ne,ni,vi,te,ti,cs,pe,pi
      CHARACTER buffer*1024

      REAL, ALLOCATABLE :: ref_pfr(:) 

      ion = 1

      output = .FALSE.

c...  Check that the reference plasma has been assigned:
      IF (ref_ntube.EQ.0) 
     .  CALL ER('InterpolateRefrencePlasma','Reference plasma does '//
     .          'not appear to be assigned',*99)

c...  Identify the tubes in the reference plasma that are on either side
c     of the focus tube (ITUBE):
      it1 = -1
      it2 = -1
      DO it = 1, ref_ntube
        IF (tube%type.NE.ref_tube(it)%type) CYCLE                  ! This may fail for more complicated double-null grids, and 
        IF (tube%rho.GE.ref_tube(it)%rho)               it1 = it   ! in particular near-connected double-null grids...
        IF (tube%rho.LT.ref_tube(it)%rho.AND.it2.EQ.-1) it2 = it
      ENDDO

      IF (output) WRITE(0,*) 'IT1,IT2',it1,it2

c...  Handle particular situations when the current focus tube is mapped
c     onto the reference grid:
      IF (it1.NE.-1.AND.it2.NE.-1) THEN
c...    All fine, do nothing:        
      ELSEIF (tube%type.EQ.GRD_SOL.AND.
     .        it1.EQ.-1.AND.it2.NE.-1) THEN
c...    In the SOL but between the separatrix and the outermost core ring:
        it1 = it2 
      ELSEIF (tube%type.EQ.GRD_CORE.AND.
     .        it1.EQ.-1.AND.it2.EQ.1) THEN
c...    Extrapolate inward from the innermost core tube:
        it1 = 2
        it2 = 1
      ELSEIF (tube%type.EQ.GRD_CORE.AND.
     .        it1.NE.-1.AND.it2.EQ.-1) THEN
c...    Between the outermost core ring and the separatrix:
        it2 = it1
      ELSEIF (tube%type.EQ.GRD_PFZ.AND.
     .        it1.NE.-1.AND.it2.EQ.-1) THEN
c...    In the PFZ but between the outermost core/PFZ tube and the separatrix:
        it2 = it1 
      ELSEIF (tube%type.EQ.GRD_PFZ.AND.
     .        it1.EQ.-1.AND.it2.NE.-1) THEN
c...    In the PFZ but beyond the innermost (farthest from separatrix) core/PFZ 
c       tube, so just take this tube as a reference ("flat" extrapolation):
        WRITE(buffer,'(A,I5)') 
     .    'Assigning plasma beyond PFZ boundary for TUBE=',tube%index
        CALL WN('InterpolateReferencePlasma',TRIM(buffer))
        it1 = it2
      ELSE
        CALL ER('InterpolateReferencePlasma','Interpolation of the '//
     .          'reference plasma failed for OSM solver option 2',*99)
      ENDIF
      
      IF (output) WRITE(0,*) 'IT1,IT2',it1,it2

      ic1 = tube%cell_index(LO)
      ic2 = tube%cell_index(HI)
      n = ic2 - ic1 + 1

      IF (output) WRITE(0,*) 'IC1,IC2',ic1,ic2,n

      ALLOCATE(ref_pfr(ref_ncell))

      DO i1 = 1, 2
c...    Interpolate along each reference ring:
        IF (i1.EQ.1) it = it1
        IF (i1.EQ.2) it = it2

        ref_ic1 = ref_tube(it)%cell_index(LO)
        ref_ic2 = ref_tube(it)%cell_index(HI)
        ref_pfr(ref_ic1:ref_ic2) = ref_cell(ref_ic1:ref_ic2)%p / 
     .                             ref_tube(it)%pmax

        IF (output) WRITE(0,*) 'REF_ICx',ref_ic1,ref_ic2
        i = 0
        DO ic = ic1, ic2
          i = i + 1
          pfr = cell(ic-ic1+1)%p / tube%pmax
c...      Identify the interpolation point on the reference ring:
          DO ref_ic = ref_ic1, ref_ic2
            int_ic1 = ref_ic - 1
            int_ic2 = ref_ic
            IF (pfr.LT.ref_pfr(ref_ic)) EXIT              
          ENDDO
          IF (ref_ic.EQ.ref_ic2+1) int_ic1 = ref_ic2

          IF (output) WRITE(0,*) 'INT_ICx',int_ic1,int_ic2
c...      Volume fluid quantities:
          IF (int_ic2.EQ.ref_ic1) THEN
c          IF (int_ic2.EQ.1) THEN  ! bug, -SL, 04/06/2013
            itarget = LO
            fr = pfr / ref_pfr(int_ic2)
            val1(1,i) = ref_tube(it)%ne(itarget)  
            val1(2,i) = ref_tube(it)%ni(itarget,ion)  
            val1(3,i) = ref_tube(it)%vi(itarget,ion)  
            val1(4,i) = ref_tube(it)%te(itarget)  
            val1(5,i) = ref_tube(it)%ti(itarget,ion)  
          ELSE
            fr = (    pfr          - ref_pfr(int_ic1)) /
     .           (ref_pfr(int_ic2) - ref_pfr(int_ic1))
            val1(1,i) = ref_fluid(int_ic1,ion)%ne 
            val1(2,i) = ref_fluid(int_ic1,ion)%ni 
            val1(3,i) = ref_fluid(int_ic1,ion)%vi 
            val1(4,i) = ref_fluid(int_ic1,ion)%te 
            val1(5,i) = ref_fluid(int_ic1,ion)%ti 
          ENDIF
          IF (int_ic1.EQ.ref_ic2) THEN
            itarget = HI
            fr = (      pfr         - ref_pfr(int_ic1)) /
     .           (ref_tube(it)%pmax - ref_pfr(int_ic1))
            val2(1,i) = ref_tube(it)%ne(itarget)  
            val2(2,i) = ref_tube(it)%ni(itarget,ion)  
            val2(3,i) = ref_tube(it)%vi(itarget,ion)  
            val2(4,i) = ref_tube(it)%te(itarget)  
            val2(5,i) = ref_tube(it)%ti(itarget,ion)  
          ELSE
            val2(1,i) = ref_fluid(int_ic2,ion)%ne 
            val2(2,i) = ref_fluid(int_ic2,ion)%ni 
            val2(3,i) = ref_fluid(int_ic2,ion)%vi 
            val2(4,i) = ref_fluid(int_ic2,ion)%te 
            val2(5,i) = ref_fluid(int_ic2,ion)%ti 
          ENDIF
          IF (i1.EQ.1) val3(1:5,i) = (1.0-fr)*val1(1:5,i)+fr*val2(1:5,i)
          IF (i1.EQ.2) val4(1:5,i) = (1.0-fr)*val1(1:5,i)+fr*val2(1:5,i)
 
c          IF (output.AND.ic.LT.ic1+5)
c     .      WRITE(0,'(A,I6,2(I6,F10.4),1P,4E10.2,0P)') 
c     .        'INT:',i1,i,pfr,int_ic1-ref_ic1+1,fr,
c     .        val1(1,i),val2(1,i),val3(1,i),val4(1,i)

          IF (output.AND.ic.LT.ic1+11)
     .      WRITE(0,'(A,I6,2(I6,F10.4),4F10.2)') 
     .        'INT:',i1,i,pfr,int_ic1-ref_ic1+1,fr,
     .        val1(5,i),val2(5,i),val3(5,i),val4(5,i)
        ENDDO
c...    Target data:
        IF (tube%type.NE.GRD_CORE) THEN
          DO itarget = LO, HI
            IF (i1.EQ.1) THEN
              val5(1,itarget) = ref_tube(it)%ne(itarget)  
              val5(2,itarget) = ref_tube(it)%ni(itarget,ion)  
              val5(3,itarget) = ref_tube(it)%vi(itarget,ion)  
              val5(4,itarget) = ref_tube(it)%te(itarget)  
              val5(5,itarget) = ref_tube(it)%ti(itarget,ion)            
              IF (output) WRITE(0,*) 'ASSIGNING VAL5',it
            ELSE
              val6(1,itarget) = ref_tube(it)%ne(itarget)  
              val6(2,itarget) = ref_tube(it)%ni(itarget,ion)  
              val6(3,itarget) = ref_tube(it)%vi(itarget,ion)  
              val6(4,itarget) = ref_tube(it)%te(itarget)  
              val6(5,itarget) = ref_tube(it)%ti(itarget,ion)            
              IF (output) WRITE(0,*) 'ASSIGNING VAL6',it
            ENDIF
          ENDDO
        ENDIF
      ENDDO

c...  Set the weight function between the interpolation tubes on 
c     either side of the focus tube:
      IF (it1.EQ.it2) THEN
        fr = 0.0
      ELSE
        fr = (    tube%rho      - ref_tube(it1)%rho) / 
     .       (ref_tube(it2)%rho - ref_tube(it1)%rho)
      ENDIF

      IF (output) THEN
        WRITE(0,*) 'IT1,IT2,FR=',it1,it2,fr
        WRITE(0,*) 'RHO(ITUBE)=',tube%rho
        WRITE(0,*) 'RHO1,2    =',ref_tube(it1)%rho,ref_tube(it2)%rho
      ENDIF
  
c...  Assign volume plasma data:
      fluid(1:n,ion)%ne = (1.0-fr) * val3(1,1:n) + fr * val4(1,1:n)
      fluid(1:n,ion)%ni = (1.0-fr) * val3(2,1:n) + fr * val4(2,1:n)
      fluid(1:n,ion)%vi = (1.0-fr) * val3(3,1:n) + fr * val4(3,1:n)
      fluid(1:n,ion)%te = (1.0-fr) * val3(4,1:n) + fr * val4(4,1:n) 
      fluid(1:n,ion)%ti = (1.0-fr) * val3(5,1:n) + fr * val4(5,1:n)
      IF (output) THEN
        DO i = 1, 10
          WRITE(0,*) 'INTER:',i,fluid(i,ion)%ne,fluid(i,ion)%te,
     .                          fluid(i,ion)%ti
        ENDDO
      ENDIF
c...  Assign target data:
      IF (tube%type.NE.GRD_CORE) THEN
        DO it = LO, HI
          ne = (1.0 - fr) * val5(1,it) + fr * val6(1,it)     
          ni = (1.0 - fr) * val5(2,it) + fr * val6(2,it)        
          vi = (1.0 - fr) * val5(3,it) + fr * val6(3,it)        
          te = (1.0 - fr) * val5(4,it) + fr * val6(4,it)     
          ti = (1.0 - fr) * val5(5,it) + fr * val6(5,it)     
          mi = 2.0 * AMU                     ! *** hardcoded: not good ***
          cs = SQRT((te + ti) * ECH / mi)    ! Needs improvement... dediated function
          pe = ne * te * ECH                 ! Same...
          pi = ne * (ti * ECH + mi * vi**2)  ! Same...  *** using electron density for now ***
c...      This list must be the same as the main list at the end
c         of SOL28_V4:
          tube%jsat       (it,ion) = cs * ne * ECH
          tube%ne         (it)     = ne
          tube%pe         (it)     = pe
          tube%te         (it)     = te
          tube%ni         (it,ion) = ni
          tube%vi         (it,ion) = vi
          tube%machno     (it)     = ABS(vi) / cs
          tube%pi         (it,ion) = pi
          tube%ti         (it,ion) = ti
          tube%gamma      (it,ion) = -1.0
          tube%qe         (it,ion) = -1.0  ! Pass from SOLPS
          tube%te_upstream(it,ion) = -1.0

          IF (output) THEN
            WRITE(0,'(A,2I6,F10.4,1P,5E10.2,0P)') 
     .        'TARGET:',it,tube%index,fr,
     .        tube%ne(it),
     .        tube%ni(it,ion),
     .        tube%vi(it,ion),
     .        tube%te(it),
     .        tube%ti(it,ion)
            WRITE(0,'(A,2I6,F10.4,1P,5E10.2,0P)') 
     .        'VAL5  :',it,tube%index,fr,val5(1:5,it)
            WRITE(0,'(A,2I6,F10.4,1P,5E10.2,0P)') 
     .        'VAL6  :',it,tube%index,fr,val6(1:5,it)
          ENDIF
        ENDDO
      ENDIF

      DEALLOCATE(ref_pfr)      

      RETURN
 99   WRITE(0,*) 'ITUBE = ',tube%index
      WRITE(0,*) 'RHO   = ',tube%rho
      WRITE(0,*) 'TYPE  = ',tube%type
      WRITE(0,*) 'IT1,2 = ',it1,it2
      STOP
      END
c
c ====================================================================
c
c subroutine: InterpolateReferencePlasma_OLD
c
c  Treats the core, SOL and PFR separately, i.e. the interpolation 
c  will make a step at the separatrix if the focus grid is narrower 
c  there than the reference grid/plasma.
c 
      SUBROUTINE InterpolateReferencePlasma_OLD(itube)
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_solps_params
      USE mod_solps
      IMPLICIT none

      INTEGER, INTENT(IN) :: itube
      

      REAL, PARAMETER :: TOL = 1.0E-03 , ECH = 1.6022E-19, 
     .                   AMU = 1.67E-27

      INTEGER ion,it,it1,it2,ref_ic,ref_ic1,ref_ic2,int_ic1,int_ic2,
     .        i1,ic,ic1,ic2,itarget,n,i
      LOGICAL output
      REAL    pfr,fr,
     .        val1(5,1000),val2(5,1000),val3(5,1000),val4(5,1000),
     .        val5(5,2),val6(5,2),
     .        mi,ne,ni,vi,te,ti,cs,pe,pi
      CHARACTER buffer*1024

      REAL, ALLOCATABLE :: ref_pfr(:) 

      ion = 1

      output = .FALSE.

      STOP 'SHOULD NOT BE HERE AT ALL'

c...  Check that the reference plasma has been assigned:
      IF (ref_ntube.EQ.0) 
     .  CALL ER('InterpolateRefrencePlasma','Reference plasma does '//
     .          'not appear to be assigned',*99)

c...  Identify the tubes in the reference plasma that are on either side
c     of the focus tube (ITUBE):
      it1 = -1
      it2 = -1
      DO it = 1, ref_ntube
        IF (tube(itube)%type.NE.ref_tube(it)%type) CYCLE                  ! This may fail for more complicated double-null grids, and 
        IF (tube(itube)%rho.GE.ref_tube(it)%rho)               it1 = it   ! in particular near-connected double-null grids...
        IF (tube(itube)%rho.LT.ref_tube(it)%rho.AND.it2.EQ.-1) it2 = it
      ENDDO

      IF (output) WRITE(0,*) 'IT1,IT2',it1,it2

c...  Handle particular situations when the current focus tube is mapped
c     onto the reference grid:
      IF (it1.NE.-1.AND.it2.NE.-1) THEN
c...    All fine, do nothing:        
      ELSEIF (tube(itube)%type.EQ.GRD_SOL.AND.
     .        it1.EQ.-1.AND.it2.NE.-1) THEN
c...    In the SOL but between the separatrix and the outermost core ring:
        it1 = it2 
      ELSEIF (tube(itube)%type.EQ.GRD_CORE.AND.
     .        it1.EQ.-1.AND.it2.EQ.1) THEN
c...    Extrapolate inward from the innermost core tube:
        it1 = 2
        it2 = 1
      ELSEIF (tube(itube)%type.EQ.GRD_CORE.AND.
     .        it1.NE.-1.AND.it2.EQ.-1) THEN
c...    Between the outermost core ring and the separatrix:
        it2 = it1
      ELSEIF (tube(itube)%type.EQ.GRD_PFZ.AND.
     .        it1.NE.-1.AND.it2.EQ.-1) THEN
c...    In the PFZ but between the outermost core/PFZ tube and the separatrix:
        it2 = it1 
      ELSEIF (tube(itube)%type.EQ.GRD_PFZ.AND.
     .        it1.EQ.-1.AND.it2.NE.-1) THEN
c...    In the PFZ but beyond the innermost (farthest from separatrix) core/PFZ 
c       tube, so just take this tube as a reference ("flat" extrapolation):
        WRITE(buffer,'(A,I5)') 
     .    'Assigning plasma beyond PFZ boundary for TUBE=',itube
        CALL WN('InterpolateReferencePlasma',TRIM(buffer))
        it1 = it2
      ELSE
        CALL ER('InterpolateReferencePlasma','Interpolation of the '//
     .          'reference plasma failed for OSM solver option 2',*99)
      ENDIF
      
      IF (output) WRITE(0,*) 'IT1,IT2',it1,it2

      ic1 = tube(itube)%cell_index(LO)
      ic2 = tube(itube)%cell_index(HI)
      n = ic2 - ic1 + 1

      ALLOCATE(ref_pfr(ref_ncell))

      DO i1 = 1, 2
c...    Interpolate along each reference ring:
        IF (i1.EQ.1) it = it1
        IF (i1.EQ.2) it = it2

        ref_ic1 = ref_tube(it)%cell_index(LO)
        ref_ic2 = ref_tube(it)%cell_index(HI)
        ref_pfr(ref_ic1:ref_ic2) = ref_cell(ref_ic1:ref_ic2)%p / 
     .                             ref_tube(it)%pmax
        i = 0
        DO ic = ic1, ic2
          i = i + 1
          pfr = cell(ic)%p / tube(itube)%pmax
c...      Identify the interpolation point on the reference ring:
          DO ref_ic = ref_ic1, ref_ic2
            int_ic1 = ref_ic - 1
            int_ic2 = ref_ic
            IF (pfr.LT.ref_pfr(ref_ic)) EXIT              
          ENDDO
          IF (ref_ic.EQ.ref_ic2+1) int_ic1 = ref_ic2
c...      Volume fluid quantities:
          IF (int_ic2.EQ.1) THEN
            itarget = LO
            fr = pfr / ref_pfr(1)
            val1(1,i) = ref_tube(it)%ne(itarget)  
            val1(2,i) = ref_tube(it)%ni(itarget,ion)  
            val1(3,i) = ref_tube(it)%vi(itarget,ion)  
            val1(4,i) = ref_tube(it)%te(itarget)  
            val1(5,i) = ref_tube(it)%ti(itarget,ion)  
          ELSE
            fr = (    pfr          - ref_pfr(int_ic1)) /
     .           (ref_pfr(int_ic2) - ref_pfr(int_ic1))
            val1(1,i) = ref_fluid(int_ic1,ion)%ne 
            val1(2,i) = ref_fluid(int_ic1,ion)%ni 
            val1(3,i) = ref_fluid(int_ic1,ion)%vi 
            val1(4,i) = ref_fluid(int_ic1,ion)%te 
            val1(5,i) = ref_fluid(int_ic1,ion)%ti 
          ENDIF
          IF (int_ic1.EQ.ref_ic2) THEN
            itarget = HI
            fr = (      pfr         - ref_pfr(int_ic1)) /
     .           (ref_tube(it)%pmax - ref_pfr(int_ic1))
            val2(1,i) = ref_tube(it)%ne(itarget)  
            val2(2,i) = ref_tube(it)%ni(itarget,ion)  
            val2(3,i) = ref_tube(it)%vi(itarget,ion)  
            val2(4,i) = ref_tube(it)%te(itarget)  
            val2(5,i) = ref_tube(it)%ti(itarget,ion)  
          ELSE
            val2(1,i) = ref_fluid(int_ic2,ion)%ne 
            val2(2,i) = ref_fluid(int_ic2,ion)%ni 
            val2(3,i) = ref_fluid(int_ic2,ion)%vi 
            val2(4,i) = ref_fluid(int_ic2,ion)%te 
            val2(5,i) = ref_fluid(int_ic2,ion)%ti 
          ENDIF
          IF (i1.EQ.1) val3(1:5,i) = (1.0-fr)*val1(1:5,i)+fr*val2(1:5,i)
          IF (i1.EQ.2) val4(1:5,i) = (1.0-fr)*val1(1:5,i)+fr*val2(1:5,i)
 
          IF (output.AND.ic.LT.ic1+5)
     .      WRITE(0,'(A,I6,2(I6,F10.4),1P,4E10.2,0P)') 
     .        'INT:',i1,i,pfr,int_ic1-ref_ic1+1,fr,
     .        val1(1,i),val2(1,i),val3(1,i),val4(1,i)
        ENDDO
c...    Target data:
        IF (tube(itube)%type.NE.GRD_CORE) THEN
          DO itarget = LO, HI
            IF (i1.EQ.1) THEN
              val5(1,itarget) = ref_tube(it)%ne(itarget)  
              val5(2,itarget) = ref_tube(it)%ni(itarget,ion)  
              val5(3,itarget) = ref_tube(it)%vi(itarget,ion)  
              val5(4,itarget) = ref_tube(it)%te(itarget)  
              val5(5,itarget) = ref_tube(it)%ti(itarget,ion)            
              IF (output) WRITE(0,*) 'ASSIGNING VAL5',it
            ELSE
              val6(1,itarget) = ref_tube(it)%ne(itarget)  
              val6(2,itarget) = ref_tube(it)%ni(itarget,ion)  
              val6(3,itarget) = ref_tube(it)%vi(itarget,ion)  
              val6(4,itarget) = ref_tube(it)%te(itarget)  
              val6(5,itarget) = ref_tube(it)%ti(itarget,ion)            
              IF (output) WRITE(0,*) 'ASSIGNING VAL6',it
            ENDIF
          ENDDO
        ENDIF
      ENDDO

c...  Set the weight function between the interpolation tubes on 
c     either side of the focus tube:
      IF (it1.EQ.it2) THEN
        fr = 0.0
      ELSE
        fr = (    tube(itube)%rho - ref_tube(it1)%rho) / 
     .       (ref_tube(it2  )%rho - ref_tube(it1)%rho)
      ENDIF

      IF (output) THEN
        WRITE(0,*) 'IT1,IT2,FR=',it1,it2,fr
        WRITE(0,*) 'RHO(ITUBE)=',tube(itube)%rho
        WRITE(0,*) 'RHO1,2    =',ref_tube(it1)%rho,ref_tube(it2)%rho
      ENDIF
  
c...  Assign volume plasma data:
      fluid(ic1:ic2,ion)%ne = (1.0-fr) * val3(1,1:n) + fr * val4(1,1:n)
      fluid(ic1:ic2,ion)%ni = (1.0-fr) * val3(2,1:n) + fr * val4(2,1:n)
      fluid(ic1:ic2,ion)%vi = (1.0-fr) * val3(3,1:n) + fr * val4(3,1:n)
      fluid(ic1:ic2,ion)%te = (1.0-fr) * val3(4,1:n) + fr * val4(4,1:n) 
      fluid(ic1:ic2,ion)%ti = (1.0-fr) * val3(5,1:n) + fr * val4(5,1:n)
c...  Assign target data:
      IF (tube(itube)%type.NE.GRD_CORE) THEN
        DO it = LO, HI
          ne = (1.0 - fr) * val5(1,it) + fr * val6(1,it)     
          ni = (1.0 - fr) * val5(2,it) + fr * val6(2,it)        
          vi = (1.0 - fr) * val5(3,it) + fr * val6(3,it)        
          te = (1.0 - fr) * val5(4,it) + fr * val6(4,it)     
          ti = (1.0 - fr) * val5(5,it) + fr * val6(5,it)     
          mi = 2.0 * AMU                     ! *** hardcoded: not good ***
          cs = SQRT((te + ti) * ECH / mi)    ! Needs improvement... dediated function
          pe = ne * te * ECH                 ! Same...
          pi = ne * (ti * ECH + mi * vi**2)  ! Same...  *** using electron density for now ***
c...      This list must be the same as the main list at the end
c         of SOL28_V4:
          tube(itube)%jsat       (it,ion) = cs * ne * ECH
          tube(itube)%ne         (it)     = ne
          tube(itube)%pe         (it)     = pe
          tube(itube)%te         (it)     = te
          tube(itube)%ni         (it,ion) = ni
          tube(itube)%vi         (it,ion) = vi
          tube(itube)%machno     (it)     = ABS(vi) / cs
          tube(itube)%pi         (it,ion) = pi
          tube(itube)%ti         (it,ion) = ti
          tube(itube)%gamma      (it,ion) = -1.0
          tube(itube)%qe         (it,ion) = -1.0  ! Pass from SOLPS
          tube(itube)%te_upstream(it,ion) = -1.0

          IF (output) THEN
            WRITE(0,'(A,2I6,F10.4,1P,5E10.2,0P)') 
     .        'TARGET:',it,itube,fr,
     .        tube(itube)%ne(it),
     .        tube(itube)%ni(it,ion),
     .        tube(itube)%vi(it,ion),
     .        tube(itube)%te(it),
     .        tube(itube)%ti(it,ion)
            WRITE(0,'(A,2I6,F10.4,1P,5E10.2,0P)') 
     .        'VAL5  :',it,itube,fr,val5(1:5,it)
            WRITE(0,'(A,2I6,F10.4,1P,5E10.2,0P)') 
     .        'VAL6  :',it,itube,fr,val6(1:5,it)
          ENDIF
        ENDDO
      ENDIF

      DEALLOCATE(ref_pfr)      

c      STOP 'Here....'

      RETURN
 99   WRITE(0,*) 'ITUBE = ',itube
      WRITE(0,*) 'RHO   = ',tube(itube)%rho
      WRITE(0,*) 'TYPE  = ',tube(itube)%type
      WRITE(0,*) 'IT1,2 = ',it1,it2
      STOP
      END
c
c ====================================================================
c
c subroutine: AssignSOLPSPlasma
c
c
c
c
      SUBROUTINE AssignSOLPSPlasma(itube)
      USE mod_sol28_params
      USE mod_sol28_global
      USE mod_solps_params
      USE mod_solps
      IMPLICIT none

      INTEGER, INTENT(IN) :: itube

      REAL, PARAMETER :: TOL = 1.0E-03 , ECH = 1.6022E-19, 
     .                   AMU = 1.67E-27

      INTEGER ion,ic1,ic2,icell,imap,solps_index(5),is,i,itarget,ic,
     .        shift
      LOGICAL mismatch_warning
      REAL    mi,ne,ni,vi,te,ti,cs,pe,pi
     
      INTEGER, ALLOCATABLE :: map_solps(:)

      DATA mismatch_warning /.FALSE./
      SAVE

      ion = 1

c...  MAP_OSM, which maps the SOLPS data index arrays to the OSM array index
c     has not been setup:
      IF (.NOT.ALLOCATED(map_osm)) THEN
        ALLOCATE(map_osm(ncell))      
        map_osm = 0
        IF (ALLOCATED(map_divimp)) THEN
c...      OSM is called from within DIVIMP, so use the index
c         mapping information that's been setup:
          DO icell = 1, ncell
            map_osm(icell) = map_divimp(cell(icell)%ik,cell(icell)%ir)
          ENDDO
        ELSE
c...      Build the map based on cell centres, on the fly (slow...):
          ALLOCATE(map_solps(solps_n))      
          map_solps = 0
          DO icell = 1, ncell 
            DO is = 1, solps_n
             IF (map_solps(is).NE.0) CYCLE
              IF (icell.EQ.1) THEN
                WRITE(logfp,'(A,2(2F14.7,2X))') 'CEN CHECK->',
     .              cell(icell)%cencar(1:2),
     .              solps_cen(is,1:2)
              ENDIF 

             IF (ABS(cell(icell)%cencar(1)-solps_cen(is,1)).LT.TOL.AND.
     .           ABS(cell(icell)%cencar(2)-solps_cen(is,2)).LT.TOL) THEN
               map_osm  (icell) = is
               map_solps(is   ) = 1
             ENDIF
            ENDDO
          ENDDO
          DEALLOCATE(map_solps)
        ENDIF
c...    Check that the map is defined for each cell:
        DO icell = 1, ncell
          IF (map_osm(icell).EQ.0) 
     .      CALL ER('AssignSOLPSPlasma','OSM_MAP assignment',*99) 
        ENDDO
      ENDIF

c...  Identify which SOLPS data sets go with what:
c
c     *** Need to do this dyamically, based on the particle species
c         that are defined for this OSM run, so that the impurities
c         calculated in SOLPS can be included here ***
c
      solps_index = 0
      DO is = 1, nsolps_data
        IF (solps_data(is)%type  .EQ.SOLPS_DENSITY    .AND.
     .      solps_data(is)%charge.EQ.-1               .AND.
     .      solps_data(is)%z     .EQ.0                .AND.
     .      solps_data(is)%a     .EQ.0            ) solps_index(1) = is
        IF (solps_data(is)%type  .EQ.SOLPS_DENSITY    .AND.
     .      solps_data(is)%charge.EQ.1                .AND.
     .      solps_data(is)%z     .EQ.1                .AND.
     .      solps_data(is)%a     .EQ.2            ) solps_index(2) = is
        IF (solps_data(is)%type  .EQ.SOLPS_VELOCITY   .AND.
     .      solps_data(is)%charge.EQ.1                .AND.
     .      solps_data(is)%z     .EQ.1                .AND.
     .      solps_data(is)%a     .EQ.2            ) solps_index(3) = is
        IF (solps_data(is)%type  .EQ.SOLPS_TEMPERATURE.AND.
     .      solps_data(is)%charge.EQ.-1               .AND.
     .      solps_data(is)%z     .EQ.0                .AND.
     .      solps_data(is)%a     .EQ.0            ) solps_index(4) = is
        IF (solps_data(is)%type  .EQ.SOLPS_TEMPERATURE.AND.
     .      solps_data(is)%charge.EQ.1                .AND.
     .      solps_data(is)%z     .EQ.1                .AND.
     .      solps_data(is)%a     .EQ.2            ) solps_index(5) = is
      ENDDO      
c...  Check that the required fluid quantities were found:
      DO i = 1, 5
        IF (solps_index(i).EQ.0) 
     .    CALL ER('AssignSOLPSPlasma','Required data not found',*99)
      ENDDO

c...  Copy SOLPS data to the OSM fluid solution arrays:
      ic1 = tube(itube)%cell_index(LO)
      ic2 = tube(itube)%cell_index(HI)
      DO icell = ic1, ic2
        imap = map_osm(icell)
        IF (ABS(cell(icell)%cencar(1)-solps_cen(imap,1)).GT.TOL.OR.
     .      ABS(cell(icell)%cencar(2)-solps_cen(imap,2)).GT.TOL) THEN
          IF (ALLOCATED(map_divimp)) THEN
            IF (.NOT.mismatch_warning) THEN
              mismatch_warning = .TRUE.
              CALL WN('AssignSOLPSPlasma','Cell position mismatch')              
            ENDIF
          ELSE
            CALL ER('AssignSOLPSPlasma','Cell position mismatch',*99)
          ENDIF
        ENDIF
        fluid(icell,ion)%ne = solps_data(solps_index(1))%data(imap)  
        fluid(icell,ion)%ni = solps_data(solps_index(2))%data(imap)  
        fluid(icell,ion)%vi = solps_data(solps_index(3))%data(imap)
        fluid(icell,ion)%te = solps_data(solps_index(4))%data(imap)
        fluid(icell,ion)%ti = solps_data(solps_index(5))%data(imap)
      ENDDO

      fluid(ic1:ic2,ion)%parion = 0.0
      fluid(ic1:ic2,ion)%parrec = 0.0
      fluid(ic1:ic2,ion)%parano = 0.0
      fluid(ic1:ic2,ion)%momsrc = 0.0
      fluid(ic1:ic2,ion)%eneano = 0.0
      fluid(ic1:ic2,ion)%eneion = 0.0

c...  Calculate the boundary conditions at the target for open 
c     field lines (SOL and PFR's):

      shift = 0
      IF (solps_indexing.EQ.1) shift = 1

      IF (tube(itube)%type.NE.GRD_CORE) THEN
        DO itarget = LO, HI
          icell = tube(itube)%cell_index(itarget)

          imap = map_osm(icell)

          IF (itarget.EQ.LO) THEN
            imap = imap - shift
          ELSE
            imap = imap + shift
          ENDIF
          mi = 2.0 * AMU  ! *** hardcoded: not good ***
          ne = solps_data(solps_index(1))%data(imap)
          ni = solps_data(solps_index(2))%data(imap)
          vi = solps_data(solps_index(3))%data(imap)
          te = solps_data(solps_index(4))%data(imap)
          ti = solps_data(solps_index(5))%data(imap)

c          ne = fluid(icell,ion)%ne
c          ni = fluid(icell,ion)%ni
c          vi = fluid(icell,ion)%vi
c          te = fluid(icell,ion)%te
c          ti = fluid(icell,ion)%ti

          cs = SQRT((te + ti) * ECH / mi)    ! Needs improvement... dediated function
          pe = ne * te * ECH                 ! Same...
          pi = ne * (ti * ECH + mi * vi**2)  ! Same...  *** using electron density for now ***

c...      This list must be the same as the main list at the end
c         of SOL28_V4:
          tube(itube)%jsat       (itarget,ion) = vi * ne * ECH ! cs * ne * ECH
          tube(itube)%ne         (itarget)     = ne
          tube(itube)%pe         (itarget)     = pe
          tube(itube)%te         (itarget)     = te
          tube(itube)%ni         (itarget,ion) = ni
          tube(itube)%vi         (itarget,ion) = vi
          tube(itube)%machno     (itarget)     = ABS(vi) / cs
          tube(itube)%pi         (itarget,ion) = pi
          tube(itube)%ti         (itarget,ion) = ti
          tube(itube)%gamma      (itarget,ion) = -1.0
          tube(itube)%qe         (itarget,ion) = -1.0  ! Pass from SOLPS
          tube(itube)%te_upstream(itarget,ion) = -1.0
        ENDDO
      ENDIF

     
      RETURN
 99   WRITE(0,*) ' SOLPS_INDEX=',solps_index
      WRITE(0,*) ' ICELL      =',icell
      WRITE(0,*) ' IMAP       =',imap
      WRITE(0,*) ' IK,IR      =',cell(icell)%ik,cell(icell)%ir
      IF (icell.GT.0.AND.imap.GT.0) THEN
        WRITE(0,*) ' CELL%CENCAR=',cell(icell)%cencar(1:2)
        WRITE(0,*) ' SOLPS_CEN  =',solps_cen(imap,1:2)
      ENDIF
      STOP
      END
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
      REAL    node_ne

      IF (nnode.NE.3.OR.mnode.NE.2) 
     .  CALL ER('AssignCorePlasma','Non-standard core node setup '//
     .          '(using the correct grid?)',*99)

      ion = 1

      ic1 = tube(itube)%cell_index(LO)
      ic2 = tube(itube)%cell_index(HI)

      IF (node(mnode)%ne.EQ.0.0) THEN
        IF (node(mnode)%pe.NE.0.0) THEN        
          node_ne = node(mnode)%pe / node(mnode)%te  ! *** flow not taken into account! ***
        ELSE
          CALL ER('AssignCorePlasma','Neither the density nor the '//
     .            'pressure are defined',*99)
        ENDIF
      ELSE
        node_ne = node(mnode)%ne
      ENDIF

      fluid(ic1:ic2,ion)%ne = node_ne
      fluid(ic1:ic2,ion)%te = node(mnode)%te
      fluid(ic1:ic2,ion)%ni = node(mnode)%ne
      fluid(ic1:ic2,ion)%ti = node(mnode)%ti(ion)

      RETURN
 99   WRITE(0,*) ' ITUBE       =',itube
      WRITE(0,*) ' NNODE,MNODE =',nnode,mnode
      STOP
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
c        WRITE(logfp,*) 'DEBUG: Calling InterpolateTeProfile'
        CALL InterpolateProfile(1)

c...    Ti:
        IF (.TRUE.) THEN
          CALL InterpolateProfile(2)
        ELSE
          ion = 1
          DO target = LO, HI
            ic1 = icbnd1(target)
            ic2 = icbnd2(target)
            IF (target.EQ.LO.AND.ic1.EQ.1    ) ic1 = 0
            IF (target.EQ.HI.AND.ic2.EQ.icmax) ic2 = icmax + 1
            SELECTCASE (opt%ti(target))
              CASE (0)
                ti(ic1:ic2,ion) = te(ic1:ic2)*DBLE(opt%ti_ratio(target))
              CASE DEFAULT
                STOP 'NO USER TI OPTION YET'
            ENDSELECT
          ENDDO
        ENDIF

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
 
      INTEGER count,ion,target,inode
      LOGICAL development
      LOGICAL, SAVE :: first_bc = .TRUE., set_cnt_bc = .FALSE.

      development = .FALSE.

      ion = 1

      IF (logop.GT.0) THEN
        WRITE(logfp,*) 
        WRITE(logfp,'(A    )') 'Local control options:'
        WRITE(logfp,'(A,L10)') 'first_bc   = ',first_bc
        WRITE(logfp,'(A,L10)') 'set_cnt_bc = ',set_cnt_bc
      ENDIF

c      node_par_mode(1:nnode) = node(1:nnode)%par_mode

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

        DO target = LO, HI
c...      Reassign options if neutrals data not available:

c          IF (opt%cflukin.EQ.2) THEN
c            opt%bc = 3
c            opt%p_ion       = 1
c            opt%p_ano       = 1
c            opt%m_mom       = 1
c            opt%m_ano       = 1
c            opt%p_ion_frac  = 100.0
c            opt%te_ano_psol = 1
c            opt%te_ano      = 1
c            opt%te_ion      = 1
c          ENDIF 

          IF (.NOT.opt%pin_data) THEN
            IF (opt%p_rec(target).EQ.2) opt%p_rec(target) = 3
            IF (opt%p_ion(target).EQ.2) THEN
              IF (opt%osm_load.EQ.0.OR.ref_nion.EQ.0) THEN
                opt%p_ion(target) = 3
                IF (opt%p_ion_frac(target).GE.0.0) 
     .            opt%p_ion_frac(target) = 0.0
              ELSE
                opt%p_ion(target) = 1
                IF (opt%p_ion_frac(target).GE.0.0) 
     .            opt%p_ion_frac(target) = 100.0
                opt_p_ion_scale(target) = 0
              ENDIF
            ENDIF
            IF (opt%te_ion(target).EQ.2) THEN
              IF (opt%osm_load.EQ.0) THEN
                opt%te_ion(target) = 3
              ELSE
                opt%te_ion(target) = 1
              ENDIF
            ENDIF
            IF (opt%ti_ion(target).EQ.2) THEN
              IF (opt%osm_load.EQ.0) THEN
                opt%ti_ion(target) = 3
              ELSE
                opt%ti_ion(target) = 1
              ENDIF
            ENDIF
c...        Sheath limited regime hack (see end of AssignNodeValues):
            IF (opt%p_ion(target).EQ.999) opt%p_ion(target) = 0
            IF (opt%m_mom(target).EQ.2  ) opt%m_mom(target) = 0
          ELSE
c...        Sheath limited regime hack (see end of AssignNodeValues):
            IF (opt%p_ion(target).EQ.999) opt%p_ion(target) = 2
          ENDIF

c...      When prescribing the ionisation source, scale to sink strength
c         if default source scaling is used (no scaling):
          IF (opt%p_ion     (target).EQ.3.AND.
     .        opt%p_ion_frac(target).EQ.100.0D0) 
     .      opt%p_ion_frac(target) = 0.0

c...      Check if sources are rescaled and target are floating:
          IF (opt%bc(target).EQ.3.AND.opt%p_ion_frac(target).NE.100.0) 
     .      CALL ER('SetupSolverOptions','Source scaling not '//       ! *** NOT SURE IF THIS IS NECESSARY ***
     .              'allowed with floating target conditions',*99)     

c...      Set whether or not scaling of the ionisation source is in effect:
c         conditions in use:
          IF (opt%bc(target).LT.0) THEN
            set_cnt_bc = .TRUE.  ! For testing upstream BC's...   *** HACK ***
            opt%bc(target) = -opt%bc(target)
          ENDIF
          opt_p_ion_scale(target) = 0
          IF (opt%bc(target).EQ.1.AND.
     .        opt%p_ion_frac(target).NE.100.0) opt_p_ion_scale(target)=1 
        ENDDO

c...    Turn of the energy transport model on the first pass, since 
c       n and v are not calculated yet (helps with stability).  Instead,
c       use a simple conduction distribution:
c        DO inode = 2, nnode-1  ! *** This didn't help...
c          IF (node_par_mode(inode).EQ.6) node_par_mode(inode) = 3
c        ENDDO

c...    Turn off energy transport model if the symmetry point temperature is too low:
        IF ((node(mnode)%te.LT.2.0.OR.
     .       node(mnode)%te.LE.node(mnode-1)%te.OR.
     .       node(mnode)%te.LE.node(mnode+1)%te).AND.
     .      node(mnode)%par_mode.EQ.6) THEN
          node(mnode)%par_mode = 3
          IF (logop.GT.0)
     .      WRITE(logfp,*) 'Symmetry point temperature too low '//
     .                     'for energy transport model, reverting '//
     .                     'to conduction transport prescription'
          WRITE(0,*) '*** TURNING OFF ENERGY MODEL ***'
        ENDIF
       
      ELSEIF (count.EQ.2) THEN
      ELSEIF (count.EQ.3) THEN  ! Should perhaps force more iterations here so that things can really settle down...
      ELSE

        IF (development.AND.set_cnt_bc.AND.count.EQ.4)
     .    cnt_boundary_conditions = .TRUE.  

        IF (cnt_boundary_conditions.AND.
     .      .NOT.cnt_super(LO).AND..NOT.cnt_super(HI)) THEN

          IF (first_bc) THEN
c            first_bc = .FALSE.
            machno(0      ,ion) = 1.0D0
            machno(icmax+1,ion) = 1.0D0
            cnt_target    = .FALSE.
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

            IF (.FALSE..AND.count.EQ.5) THEN   
              parion(icbnd1(LO):icbnd2(LO),ion) =     
     .          0.5D0 * parion(icbnd1(LO):icbnd2(LO),ion)
              eneion(icbnd1(LO):icbnd2(LO),ion) = 
     .          0.5D0 * eneion(icbnd1(LO):icbnd2(LO),ion)  

              parion(icbnd1(HI):icbnd2(HI),ion) =     
     .          1.5D0 * parion(icbnd1(HI):icbnd2(HI),ion)
              eneion(icbnd1(HI):icbnd2(HI),ion) = 
     .          1.5D0 * eneion(icbnd1(HI):icbnd2(HI),ion)  
            ENDIF

            IF (.FALSE..AND.count.EQ.4) THEN
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

        IF (opt%p_ano(LO).EQ.4.OR.opt%p_ano(HI).EQ.4) THEN
          cnt_particles = .TRUE.
          cnt_momentum  = .TRUE.
        ENDIF

c...    The target jsat is being set based on the upstream pressure, and need to keep 
c       updating the particle sources if super-sonic targets are active: 
        IF ((node(1    )%jsat(ion).LT.0.0.OR.
     .       node(nnode)%jsat(ion).LT.0.0).AND.
     .      (cnt_super(LO).OR.cnt_super(HI))) THEN
          cnt_particles = .TRUE.
        ENDIF


        IF (development.AND.count.EQ.4) THEN   ! For debugging energy profile with flow calculated
          cnt_energy       = .TRUE.
          cnt_prescription = .TRUE.
          cnt_boundary_conditions = .TRUE.
        ENDIF
c        cnt_options = .FALSE.
      ENDIF

c...  User control:
      CALL User_SetupSolverOptions(count)

c...  Output:
      IF (logop.GT.0) THEN
        WRITE(logfp,*) 
        WRITE(logfp,'(A    )') 'OSM control options:'
        WRITE(logfp,'(A,I10)') 'COUNT        = ',count
        WRITE(logfp,'(A,L10)') 'OPT%PIN_DATA = ',opt%pin_data
        WRITE(logfp,*) 
        WRITE(logfp,'(A,L10)') 'CNT_TARGET       = ',cnt_target
        WRITE(logfp,'(A,L10)') 'CNT_PARTICLES    = ',cnt_particles
        WRITE(logfp,'(A,L10)') 'CNT_MOMENTUM     = ',cnt_momentum
        WRITE(logfp,'(A,L10)') 'CNT_ENERGY       = ',cnt_energy
        WRITE(logfp,'(A,L10)') 'CNT_PRESCRIPTION = ',cnt_prescription
        WRITE(logfp,*) 
        WRITE(logfp,'(A,2I10)')   'BC          = ',opt%bc 
        WRITE(logfp,'(A,2I10)')   'SUPER       = ',opt%super
        WRITE(logfp,'(A,2I10)')   'P_REC       = ',opt%p_rec
        WRITE(logfp,'(A,2I10)')   'P_ION       = ',opt%p_ion
        WRITE(logfp,'(A,2F10.4)') 'P_ION_EXP   = ',opt%p_ion_exp
        WRITE(logfp,'(A,2I10)')   'P_ION_SCALE = ',opt_p_ion_scale
        WRITE(logfp,'(A,2F10.2)') 'P_ION_FRAC  = ',opt%p_ion_frac
        WRITE(logfp,'(A,2I10)')   'P_ANO       = ',opt%p_ano
        WRITE(logfp,'(A,2I10)')   'P_ANO_DIST  = ',opt%p_ano_dist
        WRITE(logfp,'(A,2F10.2)') 'P_ANO_EXP   = ',opt%p_ano_exp
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
        WRITE(logfp,'(A,  I10)')  'NNODE       =',nnode
c        WRITE(logfp,'(A,10I10)') 'PAR_MODE=',node_par_mode(1:nnode)
        WRITE(logfp,'(A,10I10)')  'PAR_MODE    =',node(1:nnode)%par_mode
        WRITE(logfp,*) 
      ENDIF

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE AssignTargetConditions(count)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER, INTENT(IN) :: count

      REAL*8 GetNodePressure

      INTEGER ion,target,ic1,ic2,ic,in,inode
      REAL*8  net(3),integral(10),cs,source(icmax),
     .        srcint(0:icmax),fact,p,vb

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
              ti(ic,ion) = DBLE(node(in)%ti(ion))
c              ti(ic,ion) = te(ic) * DBLE(opt%ti_ratio(target))  ! Need to check options: Ti from input file or ratio, etc...

c...          Specify target temperature everywhere on the half ring for first pass
c             through the solver, which will result in incorrect momentum sources
c             but this will be corrected on subsequent passes:
c             *** use linear interpolation to the symmetry point instead! ***
              IF (count.EQ.1) THEN
                ic1 = icbnd1(target)
                ic2 = icbnd2(target)
                te(ic1:ic2    ) = te(ic)
                ti(ic1:ic2,ion) = ti(ic,ion)
              ENDIF

              IF     (node(in)%jsat(ion).NE.0.0) THEN

                IF (node(in)%jsat(ion).LT.0.0) THEN

                  cs = DSQRT( (te(ic) + ti(ic,ion)) * ECH / mi(ion) )  ! Needs improvement... CalcCs
                  vb = cs * machno(ic,ion)
 
                  inode = INT(-node(in)%jsat(ion))
                  p     = GetNodePressure(inode,ion) / ECH
 
                  isat(ic,ion) = -(p * ECH) * vb / 
     .                 ( (te(ic) + ti(ic,ion)) * ECH + mi(ion) * vb**2 )

                  write(88,*) 'pres',inode,p
                  write(88,*) 'isat',isat(ic,ion)

c                 STOP 'well, made it here'

                ELSE
                  fact = area(ic) / ECH
                  isat(ic,ion) = -DBLE(node(in)%jsat(ion)) * fact
                ENDIF

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
              vi(ic,ion) = cs * machno(ic,ion) * tsign(target)

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

              gamma(1)   = 5.0D0
              gamma(2)   = 3.5D0
              qe(ic    ) = gamma(1) * isat(ic,ion)*ECH * te(ic    ) 
              qi(ic,ion) = gamma(2) * isat(ic,ion)*ECH * ti(ic,ion) 

c              write(0,*) 'qe,qi=',qe(ic),qi(ic,ion)

            CASE (2)
c...          Target particle and momentum from upstream cross-field 
c             and volume sources:
              STOP 'OPTION NO LONGER SUPPORTED'
              te(ic) = DBLE(node(in)%te)  ! Add some checks...
              te(ic) = DBLE(node(in)%ti(ion))
c              ti(ic,ion) = te(ic) * DBLE(opt%ti_ratio(target))  ! Not good, need to check options

              cs = DSQRT( (te(ic) + ti(ic,ion)) * ECH / mi(ion) )  ! Needs improvement... CalcCs
              vi(ic,ion) = cs * machno(ic,ion) * tsign(target)

              isat(ic,ion) =  0.0D0
              ne(ic)       = -1.0D0
              ni(ic,ion)   = -1.0D0
              pe(ic)       = -1.0D0
              pi(ic,ion)   = -1.0D0

            CASE (3)
c...          Floating jsat and Te at the target:

c             *** NO CHECK AS YET ON WHETHER THE REFERENCE SOLUTION IS DEFINED! ***

              IF (count.EQ.1) THEN
                ne(1:icmax)     = DBLE(ref_fluid(1:icmax,ion)%ne)
                ni(1:icmax,ion) = DBLE(ref_fluid(1:icmax,ion)%ni)
                vi(1:icmax,ion) = DBLE(ref_fluid(1:icmax,ion)%vi)
                te(1:icmax)     = DBLE(ref_fluid(1:icmax,ion)%te)
                ti(1:icmax,ion) = DBLE(ref_fluid(1:icmax,ion)%ti)
              ENDIF

              isat (ic,ion) = -DBLE(ref_tube%jsat (target,ion) / ECH)
              ne   (ic)     =  DBLE(ref_tube%ne   (target))
              ni   (ic,ion) =  DBLE(ref_tube%ni   (target,ion))
              vi   (ic,ion) =  DBLE(ref_tube%vi   (target,ion))
              pe   (ic)     =  DBLE(ref_tube%pe   (target))
              pi   (ic,ion) =  DBLE(ref_tube%pi   (target,ion))
              te   (ic)     =  DBLE(ref_tube%te   (target))
              ti   (ic,ion) =  DBLE(ref_tube%ti   (target,ion))
              gamma(target) =  DBLE(ref_tube%gamma(target,ion))
              qe   (ic)     =  DBLE(ref_tube%qe   (target,ion))

            CASEDEFAULT
              CALL ER('AssignTargetConditions','Unknown option',*99)
          ENDSELECT

        ENDDO  ! End of target loop
      ENDDO  ! End of ION loop

c      WRITE(logfp,*) 'ISAT1:',ion,isat(ictarg(LO),1),
c     .                            isat(ictarg(HI),1)



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
      SUBROUTINE EvaluateFluidSolution(count,MAX_ITERATIONS)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER, INTENT(IN) :: count, MAX_ITERATIONS
 
      INTEGER target,ictarget(2),ion,ic
      REAL*8  delta(2),adjust(2),
     .        h_ne(0:S28_MAXNKS+1),
     .        h_ni(0:S28_MAXNKS+1,0:S28_MAXNION),   
     .        h_vi(0:S28_MAXNKS+1,0:S28_MAXNION),   
     .        h_pe(0:S28_MAXNKS+1),                 
     .        h_pi(0:S28_MAXNKS+1,0:S28_MAXNION),   
     .        h_te(0:S28_MAXNKS+1),                 
     .        h_ti(0:S28_MAXNKS+1,0:S28_MAXNION)
      REAL*8, POINTER :: s(:)

      LOGICAL initial_run(2)

      DATA adjust /0.0D0, 0.0D0/
c      DATA delta, adjust /0.5D0, 0.5D0, 0.0D0, 0.0D0/
      SAVE

      ion = 1

      ictarget(LO) = 0
      ictarget(HI) = icmax + 1

      DO target = LO, HI

        IF (target.EQ.LO) s => sfor
        IF (target.EQ.HI) s => sbak

c        WRITE(logfp,*) 'debug: adjust start:',target,adjust(target)

c...    Check for near-target sonic transition:
        SELECTCASE (opt%super(target))
          CASE(0)
c           Do nothing:
          CASE(1)
c           Delta target Mach number and iterate:
            IF     (count.GT.MAX_ITERATIONS) THEN
c...          Reset ADJUST for next flux-tube:
              adjust = 0.0D0
              RETURN
            ELSEIF (anl_imaginary(target)) THEN
              cnt_super(target) = .TRUE. 
              IF     (count.EQ.1.OR.adjust(target).EQ.0.0D0) THEN
                IF ((s(anl_ic_super(target)).GT.0.90D0*s(icmid)).AND.
     .           ((target.EQ.LO.AND.node(1    )%jsat(ion).GE.0.0).OR.      ! Be less aggressive in turning off the super-sonic target
     .            (target.EQ.HI.AND.node(nnode)%jsat(ion).GE.0.0))) THEN   ! when jsat is based on the upstream pressure.
c                IF (s(anl_ic_super(target)).GT.0.67D0*s(icmid)) THEN
c                 Sonic transition is too close to the symmetry point,
c                 so turn off super sonic target search on the half-ring:
                  cnt_super(target) = .FALSE.                   
                  opt%super(target) = 0
                  adjust(target) =  -999.0D0
                  delta (target) =  -999.0D0
                  WRITE(logfp,*) 'TURNING OFF SUPER SONIC TARGET'
                ELSE
c                 Initialise:
                  adjust(target) =  0.5D0 ! 0.1D0
                  delta (target) =  0.1D0
                  initial_run(target) = .TRUE.
                  WRITE(logfp,*) 'INITIALIZING MACH SEARCH'
                ENDIF
              ELSEIF (adjust(target).LT.0.0D0) THEN
                adjust(target) = -0.5D0 * adjust(target)
              ELSEIF (adjust(target).GT.0.0D0.AND.
     .               initial_run(target).AND.
     .               MOD(count,4).EQ.0) THEN
c               Give things a kick:
                delta(target) = delta(target) * 2.0D0
                WRITE(logfp,*) 'KICK!'
              ENDIF
            ELSE
              initial_run(target) = .FALSE.
              IF (DABS(adjust(target)).LT.1.0D-2.OR.
     .                 delta (target) .LT.0.01D0) THEN
c...            Done, reset:
c                delta(target) = 0.5D0
                adjust(target) = 0.0D0
                cnt_super(target) = .FALSE. 
              ELSEIF (adjust(target).GT.0.0D0) THEN
                adjust(target) = -0.5D0 * adjust(target)
                cnt_super(target) = .TRUE. 
              ENDIF
c              IF     (adjust(target).GT.0.0D0) THEN
c                adjust(target) = -0.5D0 * adjust(target)
c                cnt_super(target) = .TRUE. 
c              ELSEIF (DABS(adjust(target)).LT.1.0D-3) THEN
cc...            Done, reset:
c                delta(target) = 0.5D0
c                adjust(target) = 0.0D0
c                cnt_super(target) = .FALSE. 
c                WRITE(logfp,*) 'TRING TO STOP'
c              ENDIF
            ENDIF

            IF (cnt_super(target)) THEN
              delta(target) = delta(target) * (1.0D0 + adjust(target))
              machno(ictarget(target),ion) = 1.0D0 + delta(target)
            ENDIF

            IF (logfp.GT.0)     
     .        WRITE(logfp,'(A,I4,L2,2I4,3F12.6,L2,I5)') 'SUPER: ',
     .          target,anl_imaginary(target),
     .          anl_ic_super(target),icmid,delta(target),adjust(target),
     .          machno(ictarget(target),ion),cnt_super(target),
     .          count

          CASEDEFAULT
            STOP 'NO USER ROUTINES YET'
        ENDSELECT

      ENDDO

c...  Calculated CHISQ for the plasma solution, to keep track of how things are evolving:
      IF (count.GT.1) THEN
        chisq = 0.0D0
        t_chisq = 0.0D0
        m_chisq = 0.0D0
        DO ic = 1, icmax
          chisq(1) = ((ne(ic    ) - h_ne(ic    )) / h_ne(ic    ))**2
          chisq(2) = ((ni(ic,ion) - h_ni(ic,ion)) / h_ni(ic,ion))**2
          chisq(3) = ((vi(ic,ion) - h_vi(ic,ion)) / h_vi(ic,ion))**2
          chisq(4) = ((pe(ic    ) - h_pe(ic    )) / h_pe(ic    ))**2
          chisq(5) = ((pi(ic,ion) - h_pi(ic,ion)) / h_pi(ic,ion))**2
          chisq(6) = ((te(ic    ) - h_te(ic    )) / h_te(ic    ))**2
          chisq(7) = ((ti(ic,ion) - h_ti(ic,ion)) / h_ti(ic,ion))**2

          t_chisq = t_chisq + chisq
          m_chisq = MAX(m_chisq,chisq)
 
c          t_chisq(1) =MAX(chisq(1),((ne(i    )-h_ne(i    ))/ne(i    ))**2)
c          t_chisq(2) =MAX(chisq(2),((ni(i,ion)-h_ni(i,ion))/ni(i,ion))**2)
c          t_chisq(3) =MAX(chisq(3),((vi(i,ion)-h_vi(i,ion))/vi(i,ion))**2)
c          t_chisq(4) =MAX(chisq(4),((pe(i    )-h_pe(i    ))/pe(i    ))**2)
c          t_chisq(5) =MAX(chisq(5),((pi(i,ion)-h_pi(i,ion))/pi(i,ion))**2)
c          t_chisq(6) =MAX(chisq(6),((te(i    )-h_te(i    ))/te(i    ))**2)
c          t_chisq(7) =MAX(chisq(7),((ti(i,ion)-h_ti(i,ion))/ti(i,ion))**2)
        ENDDO
        t_chisq = t_chisq / DBLE(icmax-1)
        m_chisq = m_chisq / DBLE(icmax-1)
        t_chisq(0) = SUM(t_chisq(1:7))
        m_chisq(0) = SUM(m_chisq(1:7))
c...    Do the target conditions as well...
      ENDIF
      IF (logop.GE.1) THEN
        WRITE(logfp,*) 
        WRITE(logfp,'(A,I4,1P,8E10.2,0P)') 'T_CHISQ:',count,t_chisq(0:7)
        WRITE(logfp,'(A,I4,1P,8E10.2,0P)') 'M_CHISQ:',count,m_chisq(0:7)
        WRITE(logfp,*) 
      ENDIF
c...  Store plasma solution:
      h_ne = ne
      h_ni = ni
      h_vi = vi
      h_pe = pe
      h_pi = pi
      h_te = te
      h_ti = ti

      RETURN
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE SOL28_V4(sol_option1,tube1,icmax1,cell1,pin1,fluid1,
     .                    field1,
     .                    ref_tube1,ref_nion1,ref_icmax1,ref_fluid1,  ! Do these need the "1" business?
     .                    nnode1,mnode1,node1,nion1,opt_global)
      USE mod_sol28_params
      USE mod_sol28
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER, INTENT(IN) :: sol_option1,icmax1,nnode1,mnode1,nion1,
     .                       ref_nion1,ref_icmax1
      TYPE(type_tube       ) :: tube1
      TYPE(type_cell       ) :: cell1     (icmax1)
      TYPE(type_neutral    ) :: pin1      (icmax1,nion1)
      TYPE(type_fluid      ) :: fluid1    (icmax1,nion1)
      TYPE(type_field      ) :: field1    (icmax1)
      TYPE(type_tube       ) :: ref_tube1
      TYPE(type_fluid      ) :: ref_fluid1(ref_icmax1,ref_nion1)
      TYPE(type_node       ) :: node1     (nnode1)
      TYPE(type_options_osm), INTENT(IN) :: opt_global

      INTEGER count,ion,ic
      LOGICAL cont

      INTEGER, PARAMETER :: MAX_ITERATIONS = 20

c...  Halt the code since data is almost certainly not being
c     passed properly, needs fixing but likely this is not an
c     issue for some time:
      IF (nion1.NE.1) CALL ER('SOL28_V4','NION not equal to 1',*99)

c...  Data transfer:  ...don't really like this... other way of doing it?
      sol_option = sol_option1
      nion = nion1
      tube = tube1
      icmax = icmax1  
      cell(1:icmax) = cell1(1:icmax)
      DO ion = 1, nion
        pin  (1:icmax,ion) = pin1  (1:icmax,ion)
        fluid(1:icmax,ion) = fluid1(1:icmax,ion)
      ENDDO
      IF (ref_icmax1.GT.1.AND.icmax1.EQ.ref_icmax1) THEN
        ref_tube = ref_tube1
        ref_nion = ref_nion1
c        WRITE(0,*) 'ref_nion,ref_icmax:',
c     .    ref_nion,ref_icmax1
        DO ion = 1, ref_nion
          ref_fluid(1:icmax,ion) = ref_fluid1(1:icmax,ion)
        ENDDO
      ELSE
        ref_nion = 0
      ENDIF
      IF (logop.GT.0.AND.ref_nion.EQ.0.AND.ref_nion1.NE.0) THEN
        WRITE(logfp,*)
        WRITE(logfp,*) 'Reference solution data not available '//
     .                 'for this tube because the geometry''s '
        WRITE(logfp,*) 'of the reference and current cases '//
     .                 'do not match.'
      ENDIF
      
      nnode = nnode1
      mnode = mnode1
      node(1:nnode) = node1(1:nnode) 

      opt = opt_global
      logop = opt%log
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
      
      ! jrh mod begin
      bfield(1:icmax) = field1(1:icmax)%b
      ! jrh mod end

      ai(1) = 2.0D0
      mi(1) = ai(1) * AMU
      zi(1) = 1.0D0
      ci(1) = 1.0D0

      gamma = 0.0D0

c...  Basic initialization:
      cont = .TRUE.
      count = 0

c...  Control:
      cnt_options = .TRUE.

      IF (logop.GT.0) THEN
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

      t_chisq = 1.0D0
      m_chisq = 1.0D0

c...  Main loop:
      DO WHILE (cont)
c...    Loop counter:
        count = count + 1

        IF (logop.GT.0) THEN
          WRITE(logfp,*) '==============================='
          WRITE(logfp,*) 'COUNT:',tube%ir,count
          WRITE(logfp,*) '==============================='
        ENDIF

c...    Set options for solver iteration:
c       ----------------------------------------------------------------
        IF (cnt_options     ) CALL SetupSolverOptions(count)

c...    Compile/modify ion quantities from kinetic codes:
c       ----------------------------------------------------------------
        IF (.TRUE.          ) CALL ProcessKineticIons  ! Only needs to be done once?

c...    Boundary conditions:
c       ----------------------------------------------------------------
        IF (cnt_target      ) CALL AssignTargetConditions(count)

c...    Set prescribed quantities:
c       ----------------------------------------------------------------
c        IF (cnt_energy      ) CALL AssignEnergySources 
c        IF (cnt_prescription) CALL PrescribeTemperatureProfiles  ! This goes in the solve fluid equations routine I think...

c...    Process source terms:
c       ----------------------------------------------------------------
        IF (cnt_particles   ) CALL AssignParticleSources
        IF (cnt_momentum    ) CALL AssignMomentumSources
        IF (cnt_energy      ) CALL AssignEnergySources 

c        IF (cnt_particles.OR.cnt_momentum) THEN
c          CALL IntegrateSources(2)
c        ENDIF

c...    In the case of floating boundary conditions:
c       ----------------------------------------------------------------
c        IF (cnt_boundary_conditions) THEN  ! Change to cnt_floating_target...
          CALL FloatingTargets
c          CALL UpdateTemperatureProfiles 
c        ENDIF

        IF (cnt_prescription) CALL PrescribeTemperatureProfiles  ! This goes in the solve fluid equations routine I think...

c              write(0,*) 'qe,qi B=',qe(0      ),qi(0      ,1)
c              write(0,*) 'qe,qi B=',qe(icmax+1),qi(icmax+1,1)

c...    Solve fluid equations:
c       ----------------------------------------------------------------
        CALL SolveFluidEquations

c...    Evaluate solution:
c       ----------------------------------------------------------------
        CALL EvaluateFluidSolution(count,MAX_ITERATIONS)

c...    Main iteration loop exit conditions:
c       ----------------------------------------------------------------
        IF     (count.LE.30.AND.(t_chisq(0).GT.1.0D-07.OR.
     .                           m_chisq(0).GT.1.0D-07)) THEN
          cont = .TRUE.
        ELSEIF (count.LE.MAX_ITERATIONS.AND.
     .          (opt%super(LO).NE.0.AND.cnt_super(LO).OR.
     .           opt%super(HI).NE.0.AND.cnt_super(HI))) THEN
c        ELSEIF (.FALSE..AND.(cnt_super(LO).OR.cnt_super(HI))) THEN
c...      Near-target sonic transition being processed:
          cont = .TRUE. 
        ELSEIF (count.GE.7) THEN
c          WRITE(logfp,*) 'SORRY, GIVING UP... NO CONVERGENCE'
          IF (count.GT.MAX_ITERATIONS.AND.logop.GT.1) 
     .      WRITE(logfp,*) 'MAXIMUM ITERATIONS REACHED'

          IF (logop.GE.1) 
     .      WRITE(logfp,*) 'SOLVER EXIT: COUNT.GE.7'

          cont = .FALSE.
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

          IF (logop.GE.1) 
     .      WRITE(logfp,*) 'SOLVER EXIT: DEFAULT'

          cont = .FALSE.
        ENDIF

      ENDDO ! End of main loop

      IF (logop.GE.1) 
     .  WRITE(logfp,*) 'SOLVER COUNT:',count

c...  Assign plasma quantities to tube arrays:
      ion = 1

c...  If any new items are added here they must also be 
c     added to the list at the end of AssignSOLPSPlasma:
      ic = ictarg(LO)
      tube%jsat  (LO,ion) = -SNGL(isat  (ic,ion)) * ECH  ! Change units in the solver to recover A m-2... what is MKS?
      tube%ne    (LO)     = SNGL(ne    (ic))
      tube%pe    (LO)     = SNGL(pe    (ic))
      tube%te    (LO)     = SNGL(te    (ic))
      tube%ni    (LO,ion) = SNGL(ni    (ic,ion))
      tube%vi    (LO,ion) = SNGL(vi    (ic,ion))
      tube%machno(LO)     = SNGL(machno(ic,ion))  ! Should this have an ION index or not...
      tube%pi    (LO,ion) = SNGL(pi    (ic,ion))
      tube%ti    (LO,ion) = SNGL(ti    (ic,ion))
      tube%gamma (LO,ion) = SNGL(gamma (LO))
      tube%qe    (LO,ion) = SNGL(qe    (ic))
      tube%te_upstream(LO,ion) = SNGL(te(icbnd2(LO)))
      tube%epot  (LO)     = SNGL(epot  (ic))

      ic = ictarg(HI)
      tube%jsat  (HI,ion) = -SNGL(isat  (ic,ion)) * ECH
      tube%ne    (HI)     =  SNGL(ne    (ic))
      tube%pe    (HI)     =  SNGL(pe    (ic))
      tube%te    (HI)     =  SNGL(te    (ic))
      tube%ni    (HI,ion) =  SNGL(ni    (ic,ion))
      tube%vi    (HI,ion) =  SNGL(vi    (ic,ion))
      tube%machno(HI)     =  SNGL(machno(ic,ion))  ! Should this have an ION index or not...
      tube%pi    (HI,ion) =  SNGL(pi    (ic,ion))
      tube%ti    (HI,ion) =  SNGL(ti    (ic,ion))
      tube%gamma (HI,ion) =  SNGL(gamma (HI))
      tube%qe    (HI,ion) =  SNGL(qe    (ic))
      tube%te_upstream(HI,ion) = SNGL(te(icbnd1(HI)))
      tube%epot  (HI)     =  SNGL(epot  (ic))

      DO ion = 1, nion
        fluid(1:icmax,ion)%ne = SNGL(ne(1:icmax))
        fluid(1:icmax,ion)%te = SNGL(te(1:icmax))
        fluid(1:icmax,ion)%ni = SNGL(ni(1:icmax,ion))
        fluid(1:icmax,ion)%vi = SNGL(vi(1:icmax,ion))
        fluid(1:icmax,ion)%ti = SNGL(ti(1:icmax,ion))

        fluid(1:icmax,ion)%parion =  SNGL(parion(1:icmax,ion))
        fluid(1:icmax,ion)%parrec = -SNGL(parrec(1:icmax,ion))
        fluid(1:icmax,ion)%parano =  SNGL(parano(1:icmax,ion))
        fluid(1:icmax,ion)%parusr =  SNGL(parusr(1:icmax,ion))
        fluid(1:icmax,ion)%parsrc =  SNGL(parsrc(1:icmax,ion))

        fluid(1:icmax,ion)%momvol = SNGL(momvol(1:icmax,ion))
        fluid(1:icmax,ion)%momano = SNGL(momano(1:icmax,ion))
        fluid(1:icmax,ion)%momusr = SNGL(momusr(1:icmax,ion))
        fluid(1:icmax,ion)%momsrc = SNGL(momsrc(1:icmax,ion))

        fluid(1:icmax,ion)%enerec = SNGL(enerec(1:icmax,ion))
        fluid(1:icmax,ion)%eneion = SNGL(eneion(1:icmax,ion))
        fluid(1:icmax,ion)%eneusr = SNGL(eneusr(1:icmax))
        fluid(1:icmax,ion)%eneano = SNGL(eneano(1:icmax))  ! Some inconsistency in the use of "ene"...
        fluid(1:icmax,ion)%enesrc = SNGL(enesrc(1:icmax,ion))

        fluid(1:icmax,ion)%eniion = SNGL(eniion(1:icmax,ion))
        fluid(1:icmax,ion)%eniusr = SNGL(eniusr(1:icmax,ion))
        fluid(1:icmax,ion)%eniano = SNGL(eniano(1:icmax,ion))  
        fluid(1:icmax,ion)%enisrc = SNGL(enisrc(1:icmax,ion))

c        WRITE(0,*) 'ENEION:',eneion(1:icmax,ion)

      ENDDO

c...  Pass data back through global arrays (ugly):
      tube1 = tube
      DO ion = 1, nion
        fluid1(1:icmax,ion) = fluid(1:icmax,ion)
      ENDDO

c...  Copy over the electrostatic field solver results:
      field1(1:icmax)%efield = efield(1:icmax)
      field1(1:icmax)%epot   = epot  (1:icmax)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE MainLoop(it1,it2,ikopt)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      INTEGER it1,it2,ikopt,cnt

      INTEGER itube,cind1,cind2,ref_cind1,ref_cind2,ref_itube,iopt,
     .        sol_option,isol,ipass,iloop
      LOGICAL cont

      LOGICAL CheckIndex

      INTEGER ion,i1,i2,i,ic1,ic2,itarget,opt_iopt
c      REAL    totsrc
      INTEGER                 nlist
      INTEGER, ALLOCATABLE :: ilist(:)

      TYPE(type_tube )              :: tube_tmp
      TYPE(type_fluid), ALLOCATABLE :: fluid_tmp(:,:)

      INTEGER nnode,mnode
      TYPE(type_node) :: node(S28_MAXNNODE)
      TYPE(type_options_osm) :: opt_tube

      INTEGER debug_cnt
      DATA    debug_cnt /0/
      SAVE

      debug_cnt = debug_cnt + 1

      IF (logop.GT.0) THEN
        WRITE(logfp,*)
        WRITE(logfp,*) 'RANGE:',it1,it2,grid%isep
        WRITE(logfp,*) 'TE   :',tube(it1)%te(LO),tube(it1)%te(HI)
      ENDIF

c...  Build list of tubes, and put the PFZ tubes at the end, and in reverse order.  This
c     helps a lot when tubes in the PFZ reference tubes that are closer to the separatrix
c     but have a higher tube index:
      ALLOCATE(ilist(it2-it1+1))
      nlist = 0
      DO ipass = 1, 2
        IF (ipass.EQ.1) THEN 
          DO itube = it1, it2
            IF (tube(itube)%type.NE.GRD_PFZ) THEN
              nlist = nlist + 1
              ilist(nlist) = itube
            ENDIF  
          ENDDO
        ELSE
          DO itube = it2, it1, -1
            IF (tube(itube)%type.EQ.GRD_PFZ) THEN
              nlist = nlist + 1
              ilist(nlist) = itube
            ENDIF  
          ENDDO
        ENDIF
      ENDDO

      ion = 1

      IF (.TRUE.) CALL ListTargetData(logfp,'Before calling solver')

      tube2(it1:it2)%state = ibclr(tube2(it1:it2)%state,1)  ! Flag that solution has not been calculated yet

c      write(0,*) 'state',it1,it2
c      write(0,*) tube2(1:it2)%state

      opt%cosm = 0

      cnt = 0
      cont  = .TRUE.
      DO WHILE (cont) 
        cnt = cnt + 1
        cont  = .FALSE.

c...    Calcualte derived quantities (efield, drifts, etc.):
c        CALL CalculateEfield
c        CALL CalculateGradients
c        CALL CalculateDrifts

c        DO itube = it1, it2
        DO iloop = 1, it2-it1+1

          itube = ilist(iloop)

          IF (ibits(tube2(itube)%state,0,1).EQ.0.AND.     ! Default symmetry point was not applied
     .        ibits(tube2(itube)%state,1,1).EQ.1) CYCLE   ! Solution has been calculated already

c...      Setup solver options:
c         CALL SetupLocalOptions(itube,opt)

          opt_iopt = 1
          IF (.TRUE.) THEN 
            DO iopt = 1, nopt
c              WRITE(0,*) 'CHECK:',opt%cflukin,
c     .                   opt_iteration(iopt)%iteration(1:2)
c              WRITE(0,*) '     :',itube,
c     .                   opt_iteration(iopt)%tube    (1:2)

c              write(0,*) 'checking --- ',itube,
c     .             CheckIndex(itube,0,opt_iteration(iopt)%tube),
c     .             opt_iteration(iopt)%tube

              IF (opt%cflukin.GE.opt_iteration(iopt)%iteration(1).AND.
     .            opt%cflukin.LE.opt_iteration(iopt)%iteration(2).AND.
     .            CheckIndex(itube,0,opt_iteration(iopt)%tube)) THEN
c     .            itube      .GE.opt_iteration(iopt)%tube(1)     .AND.
c     .            itube      .LE.opt_iteration(iopt)%tube(2)) THEN
                opt_tube = opt_iteration(iopt)
                opt_iopt = iopt
c                WRITE(0,*) 'SELECTING OPTION SET   :',itube,opt%cflukin,
c     .                     iopt
              ENDIF
            ENDDO
          ELSE
            opt_tube = opt
          ENDIF

c          WRITE(0,*) 'CHECKING SOLVER OPTION:',opt_tube%sol_n
          sol_option = 0
          DO isol = 1, opt_tube%sol_n
c            WRITE(0,*) 'CHECKING SOLVER OPTION:',isol,itube,
c     .                 opt_tube%sol_tube(1:2,isol)
            IF (itube.GE.opt_tube%sol_tube(1,isol).AND.
     .          itube.LE.opt_tube%sol_tube(2,isol)) THEN
              sol_option = opt_tube%sol_option(isol)
c              WRITE(0,*) 'SELECTING SOLVER OPTION:',isol,sol_option
            ENDIF
          ENDDO

          IF (logop.GT.0) THEN
            WRITE(logfp,'(A   )') '------------------------------------'
            WRITE(logfp,'(A,I6)') ' SOLVING TUBE: ',itube
c            WRITE(0    ,'(A,I6)') ' SOLVING TUBE: ',itube
            WRITE(logfp,'(A,I6)') ' SOL OPTON   : ',sol_option
            WRITE(logfp,'(A,I6)') ' CFLUKIN     : ',opt%cflukin
            WRITE(logfp,'(A,I6)') ' IOPT        : ',opt_iopt       
            WRITE(logfp,'(A,I6)') ' COSM        : ',opt%cosm
            WRITE(logfp,'(A,I6)') ' DEBUG_CNT   : ',debug_cnt
            WRITE(logfp,'(A   )') '------------------------------------'
          ENDIF

c...      Identify cell based data to transfer to solver:
          cind1 = tube(itube)%cell_index(LO)
          cind2 = tube(itube)%cell_index(HI)

          ref_itube = MIN(ref_ntube ,itube)
          ref_cind1 = MIN(ref_nfluid,cind1)
          ref_cind2 = MIN(ref_nfluid,cind2)
 
c...      Calculate plasma solution:
          IF (.FALSE.) THEN
c...        MPI:
          ELSE
            store_sopt (itube) = sol_option  ! *** TEMP ***

            SELECTCASE (sol_option) 
c             ----------------------------------------------------------
              CASE(1)
c...            Assign reference solution:
                IF (cind1.NE.ref_cind1.OR.cind2.NE.ref_cind2.OR.
     .              itube.NE.ref_itube) 
     .            CALL ER('MainLoop','Reference solution index '//
     .                    'mismatch',*99)
                tube(itube) = ref_tube(itube)
                DO i1 = 1, nion
                  fluid(cind1:cind2,i1) = ref_fluid(cind1:cind2,i1)
                ENDDO
c             ----------------------------------------------------------
              CASE(2)
c...            Assign SOLPS solution:
                CALL AssignSOLPSPlasma(itube)
                tube2(itube)%state = ibset(tube2(itube)%state,1)            ! Flag that solution for ITUBE has been calculated
c             ----------------------------------------------------------
              CASE(3)
c...            Interpolate reference solution:
                CALL InterpolateReferencePlasma(tube(itube),nion,
     .                 fluid(cind1:cind2,1:nion),
     .                 cell (cind1:cind2))
c                CALL InterpolateReferencePlasma(itube)
                tube2(itube)%state = ibset(tube2(itube)%state,1)            ! Flag that solution for ITUBE has been calculated
c             ----------------------------------------------------------
              CASE(28:30)
c...            SOL28 - SimpleAsPie analytic particle and momentum solver (SL)
c                  29 - CestDuGateau, the Runge-Kutta solver (JRH)
c                  30 - the time-dependent solver (JRH)
                CALL SetTargetConditions(itube)

c...            Assign solution parameter nodes:
                IF     (opt%s28mode.EQ.5.0) THEN 
                  CALL osmSetNodeValues(itube,nnode,mnode,node,opt_tube)
                ELSEIF (opt%s28mode.EQ.4.1) THEN 
                  CALL AssignNodeValues_New(itube,nnode,mnode,node,
     .                                      opt_tube)
                  
c                  WRITE(0,*) '_state:',tube2(itube)%state,
c     .                 ibits(tube2(itube)%state,0,1),
c     .                 ibits(tube2(itube)%state,1,1)

                ELSE
                  STOP 'LEGACY NODES NO LONGER SUPPORTED'
c                  CALL AssignNodeValues_Legacy(itube,nnode,mnode,node)
                ENDIF

                IF (tube(itube)%type.EQ.GRD_CORE) THEN
                  CALL AssignCorePlasma(itube,nnode,mnode,node)
                ELSE
                  CALL SOL28_V4(sol_option,tube(itube),tube(itube)%n,     ! Likely need to pass index range for 
     .                          cell (cind1:cind2),                       ! each array individually in case
     .                          pin  (cind1:cind2,nion),                  ! some have only nominal allocations, i.e.
     .                          fluid(cind1:cind2,nion),                  ! they are not in use...
     .                          field(cind1:cind2),
     .                          ref_tube(ref_itube),
     .                          ref_nion,ref_cind2-ref_cind1+1,
     .                          ref_fluid(ref_cind1:ref_cind2,ref_nion),  ! Clumsy...
     .                          nnode,mnode,node,nion,opt_tube)           ! Also: pass local options
     .                          

                  store_nnode(itube) = nnode  ! *** TEMP ***
                  store_mnode(itube) = mnode
                  store_node (1:nnode,itube) = node(1:nnode)
                  IF (store_node(1    ,itube)%jsat(ion).NE.  ! Need to save the correct jsat in case a code/flag for dynamic 
     .                tube(itube)%jsat(LO,ion))              ! assignment was used
     .                store_node(1    ,itube)%jsat(ion) =
     .                tube(itube)%jsat(LO,ion)
                  IF (store_node(nnode,itube)%jsat(ion).NE.
     .                tube(itube)%jsat(HI,ion)) 
     .                store_node(nnode,itube)%jsat(ion) =
     .                tube(itube)%jsat(HI,ion)


                ENDIF
                tube2(itube)%state = ibset(tube2(itube)%state,1)            ! Flag that solution for ITUBE has been calculated
c             ----------------------------------------------------------
              CASE DEFAULT
                WRITE(0,*) 'SOL_OPTION',sol_option
                CALL ER('MainLoop','Solver option not identified',*99)
            ENDSELECT


c...        Check to see if a portion of the OSM solution is to be over-written 
c           by the reference plasma solution:

            DO i = 1, nnode
              IF (node(i)%par_mode.NE.7) CYCLE
                
              IF (ref_ntube.EQ.0.OR.ref_ntube.EQ.1) 
     .         CALL ER('MainLoop','Reference solution overwrite found'//
     .                 ', but no reference available',*99)     

              ic1 = tube(itube)%cell_index(LO)
              ic2 = tube(itube)%cell_index(HI)     
              ALLOCATE(fluid_tmp(ic2-ic1+1,nion))
              tube_tmp = tube(itube)
              CALL InterpolateReferencePlasma(tube_tmp,nion,
     .               fluid_tmp(1:ic2-ic1+1,1:nion),cell(ic1:ic2))         
              EXIT
            ENDDO

            IF (ALLOCATED(fluid_tmp)) THEN
c              write(0,*) 'going for it!'

              DO i = 2, mnode
                IF (node(i)%par_mode.NE.7) CYCLE

c              write(0,*) 'going for it 1!',i,mnode

 
 

                i1 = node(i-1)%icell + 1
                i2 = node(i  )%icell - 1


c                write(0,*) 'numbers',ic1,ic2,i1,i2

                fluid(ic1+i1-1:ic1-1+i2,ion)%ne=fluid_tmp(i1:i2,ion)%ne
                fluid(ic1+i1-1:ic1-1+i2,ion)%vi=fluid_tmp(i1:i2,ion)%vi
                fluid(ic1+i1-1:ic1-1+i2,ion)%te=fluid_tmp(i1:i2,ion)%te
                fluid(ic1+i1-1:ic1-1+i2,ion)%ti=fluid_tmp(i1:i2,ion)%ti
              ENDDO
              DO i = mnode, nnode-1
                IF (node(i)%par_mode.NE.7) CYCLE

c              write(0,*) 'going for it 2!',i,mnode

                i1 = node(i  )%icell + 1 
                i2 = node(i+1)%icell - 1

c                write(0,*) 'numbers',ic1,ic2,i1,i2

                fluid(ic1+i1-1:ic1+i2-1,ion)%ne=fluid_tmp(i1:i2,ion)%ne
                fluid(ic1+i1-1:ic1+i2-1,ion)%vi=fluid_tmp(i1:i2,ion)%vi
                fluid(ic1+i1-1:ic1+i2-1,ion)%te=fluid_tmp(i1:i2,ion)%te
                fluid(ic1+i1-1:ic1+i2-1,ion)%ti=fluid_tmp(i1:i2,ion)%ti
              ENDDO
              DEALLOCATE(fluid_tmp)
            ENDIF


          ENDIF  ! MPI selection
        ENDDO  ! Tube loop

c        WRITE(0,*) '_count',
c     .             COUNT(ibits(tube2(it1:it2)%state,0,1).EQ.0)
c        WRITE(0,*) ibits(tube2(it1:it2)%state,0,1)

c          write(0,*) 'tube,cnt,2.1',it1,it2,cnt,
c     .                COUNT(osmnode(2:osmnnode)%type.EQ.2.1)
c          DO  i1 = it1, it2
c            write(0,*) '   state',i1,ibits(tube2(i1)%state,2,1),
c     .                               ibits(tube2(i1)%state,0,1)  
c          ENDDO


        IF    (cnt.GE.3) THEN

        ELSEIF (COUNT(osmnode(2:osmnnode)%type       .EQ.2.1).GT.0.AND.
     .          COUNT(ibits(tube2(it1:it2)%state,0,1).EQ.1  ).GT.0) THEN
c...      Check if default symmetry point was specified:
          cont = .TRUE.
        ELSEIF (COUNT(ibits(tube2(it1:it2)%state,2,1).EQ.1  ).GT.0) THEN  
c...      Check if any invalid links were present:
          IF (cnt.EQ.2) STOP 'Problemo man'
c          tube2(it1:it2)%state = IBCLR(tube2(it1:it2)%state,2)
          cont = .TRUE.
        ENDIF

c...    Modify conditions for iteration:
        CALL User_MainLoop(cont,cnt)

        tube2(it1:it2)%state = IBCLR(tube2(it1:it2)%state,2)  ! Clear registered links to undefined tubes

      ENDDO ! Iteration loop

c...  Store solution for subsequent call to solver:   *** DON'T ALWAYS NEED THIS ***
c      CALL AssignReferenceSolution  ! TURNED OFF WHILE WORKING ON THE PLASMA INTERPOLATION CODE
c                                      THIS WAS BRUTAL ANYWAY, SO NEED TO FIND A NEW SOLUTION FOR IDENTIFYING WHEN
c                                      THE REFERENCE SOLUTION SHOULD BE STORED... NOTE THAT SOME OLD CASES WILL BE 
c                                      BROKEN, BUT MOST ARE IN THE DEVELOPMENT STAGES ANYWAY...

c      WRITE(logfp,*) 
c      ion = 1
c      totsrc = 0.0
c      DO i1 = 1, ncell
c        totsrc = totsrc + fluid(i1,ion)%parsrc * cell(i1)%vol
c      ENDDO
c      WRITE(logfp,*) 'TOTSRC 3=',totsrc
c      WRITE(logfp,*)

      CALL ListTargetData(logfp,'After calling solver')

      CALL DumpData_OSM('output.solver','Done calling fluid solver')

c      CALL OutputData(logfp,'FLUID SOLUTION COMPLETE')

      DEALLOCATE(ilist)

      RETURN
 99   WRITE(0,*) '  ITUBE=',itube
      STOP
      END
