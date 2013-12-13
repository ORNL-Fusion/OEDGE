c     -*Former Mode Specification*-
c
c ======================================================================
c
      SUBROUTINE SolveFluidEquations
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER ic,ion
      REAL*8  cs

      ion = 1

      cs = -999.0  ! -CU compiler flag

      IF (nion.NE.1)
     .  CALL ER('SolveFluidEquations','One fluid ion only please',*99)

      IF (logop.GT.0) THEN
        WRITE(logfp,*) 'ABOUT TO CALL PLASMA SOLVER'
        WRITE(logfp,*) 
        WRITE(logfp,*) 'ION DATA :',mi(ion),ai(ion)
        WRITE(logfp,*) '         :',zi(ion),ci(ion)
        WRITE(logfp,*) 
        ic = 0
        WRITE(logfp,*) 'TARGET LO:',ic
        WRITE(logfp,*) '  cs,M      :',cs,machno(ic,ion)
        WRITE(logfp,*) '  jsat,area :',isat(ictarg(LO),ion),area(ic)
        WRITE(logfp,*) '  ne,ni     :',ne(ic),ni(ic,ion)
        WRITE(logfp,*) '  vi        :',vi(ic,ion)
        WRITE(logfp,*) '  pe,pi     :',pe(ic),pi(ic,ion)
        WRITE(logfp,*) '  te,ti     :',te(ic),ti(ic,ion)
        WRITE(logfp,*) '  qe,qi     :',qe(ic),qi(ic,ion)
        ic = icmax + 1
        WRITE(logfp,*) 'TARGET HI:',ic
        WRITE(logfp,*) '  cs,M      :',cs,machno(ic,ion)
        WRITE(logfp,*) '  jsat,area :',isat(ictarg(HI),ion),area(ic)
        WRITE(logfp,*) '  ne,ni     :',ne(ic),ni(ic,ion)
        WRITE(logfp,*) '  vi        :',vi(ic,ion)
        WRITE(logfp,*) '  pe,pi     :',pe(ic),pi(ic,ion)
        WRITE(logfp,*) '  te,ti     :',te(ic),ti(ic,ion)
        WRITE(logfp,*) '  qe,qi     :',qe(ic),qi(ic,ion)

c        WRITE(logfp,*) 'MOMINT:',momint(1:icmax,ion)
      ENDIF

c...  Assign target conditions to standard fluid arrays:
      

      DO ic = 1, icmax
        par(ic,ion) = parsrc(0,ion) + parint(ic,ion)
        mom(ic,ion) = momsrc(0,ion) + momint(ic,ion)
      ENDDO


      SELECTCASE (sol_option)
        CASE (28)
          CALL SimpleAsPie(ion)
c       JRH MOD BEGIN
        CASE (29)
          CALL SOLF1D_Interface(ion)
c       JRH MOD END
        
c        CASE (1)
c...      Portal to SOL23:
c         CALL BranchToSOL23
        CASEDEFAULT
          STOP 'NO CUSTOM SOLVERS YET'
      ENDSELECT


c...  For now:
      ne(1:icmax) = ni(1:icmax,ion)
      
c     jrhmod begin
c     TODO: Add input file option to control this?
      CALL Solve_Potential(ion)
c     jrhmod end


      RETURN
 99   WRITE(0,*) '  NION=',nion
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE SimpleAsPie(ion)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

      INTEGER ion

      INTEGER ic
      LOGICAL cont
      REAL*8  par1x,mom1x,a1,b1,c1,r1,cs


      cont = .TRUE.

      IF (logop.GT.0) WRITE(logfp,*) 'ICMAX:',icmax


      anl_imaginary = .FALSE.

      IF (logop.GT.0) THEN
        WRITE(logfp,*) 'SOL28'
        WRITE(logfp,'(3A4,2A10,A6,3A10,5X,2(2A10,2X),4A10)')
     .    'i','iSL','iSH','ni','vi','M','pe+pi','te','ti',
     .    'par','mom','qcond','qconv','a','b','c','r'
      ENDIF

c...  Density and velocity:
      DO WHILE (cont)
        cont = .FALSE.
        IF (logop.GT.0) THEN
          ic = 0
          WRITE(logfp,'(I4,8X,1P,D10.2,D10.2,0P,F6.2,1P,D10.2,
     .                  0P,2F10.4,5X,1P,2D10.2,2X,2D10.2)')
     .      ic,
     .      ni(ic,ion),vi(ic,ion),machno(ic,ion),pe(ic)+pi(ic,ion),
     .      te(ic),ti(ic,ion),
     .      parsrc(ic,ion),momsrc(ic,ion),
     .      qcond(ic),qconv(ic)
c          WRITE(logfp,'(A,I4,8X,42X,1P,2D10.2,2X,5D10.2,0P,3F8.2)')
c     .      '  SOL28:',ic,
c     .      parsrc(ic,ion),momsrc(ic,ion),
c     .      ni(ic,ion),vi(ic,ion),pe(ic)+pi(ic,ion),
c     .      qcond(ic),qconv(ic),te(ic),ti(ic,ion),
c     .      machno(ic,ion)
        ENDIF

        DO ic = 1, icmax

          par1x = par(ic,ion)
          mom1x = mom(ic,ion)

          a1 = (te(ic) + ti(ic,ion)) * ECH
          b1 = -mom1x
          c1 = AMU * 2.0D0 * par1x**2
          r1 = b1**2 - 4.0D0 * a1 * c1
c          b1 = -(p0 + mom1x) 
c          c1 = AMU * crmb * par1x**2

          IF (r1.GE.0.0D0) THEN
c...        Real solution:
            IF ((opt%super(LO).NE.0.AND.ic.LT.anl_ic_super(LO)).OR.
     .          (opt%super(HI).NE.0.AND.ic.GT.anl_ic_super(HI))) THEN
              ni(ic,ion) = (-b1 - DSQRT(r1)) / (2.0D0 * a1)  ! Supersonic
            ELSE
              ni(ic,ion) = (-b1 + DSQRT(r1)) / (2.0D0 * a1)
            ENDIF
          ELSE
c...        Imaginary, approximate density:
            ni(ic,ion) = -b1 / (2.0D0 * a1)
c...        Record extent of super-sonic regions:
            IF     (ic.LE.icmid) THEN
c             Only acknowledge the super-sonic flow if it is towards the 
c             appropriate target:
              IF (ic.EQ.1.OR.
     .            (c1.GE.2.AND.vi(MAX(1,ic-1),ion).LT.0.0D0)) THEN
                anl_ic_super (LO) = ic
                anl_imaginary(LO) = .TRUE.
              ENDIF
            ELSEIF (ic.GT.icmid.AND..NOT.anl_imaginary(HI)) THEN
              IF (vi(ic-1,ion).GT.0.0D0) THEN 
                anl_ic_super (HI) = ic
                anl_imaginary(HI) = .TRUE.
              ENDIF
            ENDIF
          ENDIF

c...      Set fluid velocity:
          vi(ic,ion) = par1x / ni(ic,ion)

c...      Set fluid pressure:
          ne(ic)     = ni(ic,ion)
          pe(ic)     = ne(ic) * te(ic) * ECH
          pi(ic,ion) = ni(ic,ion) * 
     .                 (ti(ic,ion) * ECH + mi(ion) * vi(ic,ion)**2)

c...      
          cs = DSQRT( (te(ic) + ti(ic,ion)) * ECH / mi(ion) )  ! Need improved calculation...
          machno(ic,ion) = DABS(vi(ic,ion)) / cs

          IF (logop.GT.0) 
     .       WRITE(logfp,'(3I4,1P,2D10.2,0P,F6.2,1P,D10.2,
     .                     0P,2F10.4,5X,1P,2D10.2,2X,2D10.2,2X,4D10.2)')
     .        ic,anl_ic_super(LO),anl_ic_super(HI),
     .        ni(ic,ion),
     .        vi(ic,ion),machno(ic,ion)*DSIGN(1.0D0,vi(ic,ion)),
     .        pe(ic)+pi(ic,ion),
     .        te(ic),ti(ic,ion),
     .        par(ic,ion),mom(ic,ion),
     .        qcond(ic),qconv(ic),
     .        a1,b1,c1,r1
c     .       WRITE(logfp,'(A,3I4,1P,4D10.2,2X,2D10.2,2X,5D10.2,
c     .                     0P,3F8.2)')
c     .        '  SOL28:',ic,anl_ic_super(LO),anl_ic_super(HI),
c     .        a1,b1,c1,r1,
c     .        par(ic,ion),mom(ic,ion),
c     .        ni(ic,ion),vi(ic,ion),pe(ic)+pi(ic,ion),
c     .        qcond(ic),qconv(ic),te(ic),ti(ic,ion),
c     .        machno(ic,ion)

        ENDDO  ! End of IC loop

        IF (logop.GT.0) THEN
          ic = icmax + 1
          WRITE(logfp,'(I4,8X,1P,D10.2,D10.2,0P,F6.2,1P,D10.2,
     .                  0P,2F10.4,5X,1P,2D10.2,2X,2D10.2)')
     .      ic,
     .      ni(ic,ion),vi(ic,ion),machno(ic,ion),pe(ic)+pi(ic,ion),
     .      te(ic),ti(ic,ion),
     .      parsrc(ic,ion),momsrc(ic,ion),
     .      qcond(ic),qconv(ic)
c          WRITE(logfp,'(A,I4,8X,42X,1P,2D10.2,2X,5D10.2,0P,3F8.2)')
c     .      '  SOL28:',ic,
c     .      parsrc(ic,ion),momsrc(ic,ion),
c     .      ni(ic,ion),vi(ic,ion),pe(ic)+pi(ic,ion),
c     .      qcond(ic),qconv(ic),te(ic),ti(ic,ion),
c     .      machno(ic,ion)
        ENDIF

      ENDDO  ! End of main loop



 
      RETURN
 99   STOP
      END
c
c ======================================================================
c
c    jrh mod begin
      SUBROUTINE SOLF1D_Interface(ion)
      USE mod_sol28_params
      USE mod_sol28_solver
      USE SOLF1D

      IMPLICIT NONE

      ! Boundary conditions
      REAL*8 BC_Dirichlet(8)

      !sources
      REAL*8 Sext_n(icmax+2),Sext_Te(icmax+2),
     .       Sext_Ti(icmax+2),Sext_v(icmax+2)

      REAL*8 Seir_n(icmax+2),Seir_Te(icmax+2),
     .       Seir_Ti(icmax+2),Seir_v(icmax+2)

      !grid
      REAL*8 x_n(icmax+2),x_v(icmax+2)
 
!       !magnetic field
      REAL*8 Bsol(icmax+2)

      !solved quantities
      REAL*8 f_n(icmax+2),f_v(icmax+1),f_Te(icmax+2),
     .       f_Ti(icmax+2),f_n0(icmax+2),f_v0(icmax+1)

      REAL*8 R, alphae, alphai, beta

      INTEGER ion, n, IT, ic, A, NEUT, IT_GRAPH, esrc

      R = 0.975D0         ! Recycling coefficient (not used)
      alphae=10000.0D0    ! Heat flux limiter
      alphai=10000.0D0    ! Heat flux limiter
      beta=0.5D0          ! Viscous flux limiter
      A = 2               ! Atomic mass number
      NEUT = 0            !NEUT=1 model with neutrals, NEUT=0 model without neutrals, if NEUT=0, recycling
                                                        !coefficient will be set to 0
      IT_GRAPH = 0        ! Do not use external graphics
      esrc     = 0        ! Solve energy equations (1) or not (0)

c    Number of grid cells
      n = icmax+1

c    Number of time steps (this should be replaced by a proper convergence check!)
      IT = 100000

c    Set up the boundary conditions
c    Density
      BC_Dirichlet(1) = ne(0)
      BC_Dirichlet(2) = ne(icmax+1)
c    Electron temperature
      BC_Dirichlet(3) = te(0)
      BC_Dirichlet(4) = te(icmax+1)
c    Ion temperature
      BC_Dirichlet(5) = ti(0,ion)
      BC_Dirichlet(6) = ti(icmax+1,ion)
c    Ion velocity
      BC_Dirichlet(7) = vi(0,ion)
      BC_Dirichlet(8) = vi(icmax+1,ion)

c    ONLY USE THIS IF ENERGY EQUATIONS ARE NOT BEING SOLVED
      IF (esrc .EQ. 0) THEN
        f_Te(1)=BC_Dirichlet(3)
        f_Te(n+1)=BC_Dirichlet(4)

        DO ic = 2, n
          f_Te(ic) = te(ic-1)
          f_Ti(ic) = ti(ic-1,ion)
        ENDDO

        f_Ti(1)=BC_Dirichlet(5)
        f_Ti(n+1)=BC_Dirichlet(6)
      ENDIF

c    Set up the sources
      DO ic=1,n+1
        Sext_n(ic)  = 0.0
        Sext_Te(ic) = 0.0
        Sext_Ti(ic) = 0.0
        Sext_v(ic)  = 0.0

        Seir_n(ic)  = 0.0
        Seir_Te(ic) = 0.0
        Seir_Ti(ic) = 0.0
        Seir_v(ic)  = 0.0
      ENDDO

c    Particle sources
      DO ic=2,n		
        Sext_n(ic) = parsrc(ic-1,ion)
      END DO
      Sext_n(1)=0.0D0
      Sext_n(n+1)=0.0D0
c    Momentum sources
      DO ic=2,n		
        Sext_v(ic) = momsrc(ic-1,ion)
      END DO
      Sext_v(1)=0.0D0
      Sext_v(n+1)=0.0D0
      
      IF (esrc .NE. 0) THEN
c       Electron energy source
        DO ic=2,n		
          Sext_Te(ic) = enesrc(ic-1,ion)
        END DO
        Sext_Te(1)=0.0D0
        Sext_Te(n+1)=0.0D0
      
c       Ion energy source
        DO ic=2,n		
          Sext_Ti(ic) = enisrc(ic-1,ion)
        END DO
        Sext_Ti(1)=0.0D0
        Sext_Ti(n+1)=0.0D0
      ENDIF

c    Set up the magnetic field
      DO ic=2,n		
        Bsol(ic) = bfield(ic-1) 
      END DO
      ! The following two lines are a fudge
      Bsol(1)=Bsol(2)
      Bsol(n+1)=Bsol(n)

c    Set up the grid 
      DO ic=1,n-1		
        x_v(ic) = cell(ic)%sbnd(1)
      END DO 	
		
      x_v(n) = cell(icmax)%sbnd(2) 	
	
      x_n(1)=x_v(1)
      x_n(n+1)=x_v(n)
      DO ic=2,n
        x_n(ic)=0.5D0*(x_v(ic-1)+x_v(ic))
      END DO
      
c      WRITE(0,*) 'IN SOLF1D, tube:', tube%ir

      CALL STEADY_STATE(f_n,f_v,f_Te,f_Ti,f_n0,f_v0,x_n,x_v,n,
     .                  Sext_n,Sext_v,Sext_Te,Sext_Ti,
     .                  Seir_n,Seir_v,Seir_Te,Seir_Ti,
     .                  Bsol,alphae,
     .                  alphai,beta,NEUT,R,A,IT,IT_graph,
     .                  BC_Dirichlet,esrc)

      ! Write the information into OSM variables
      DO ic = 1, icmax
        ni(ic,ion) = f_n(ic+1)
        vi(ic,ion) = f_v(ic+1)
        ! comment this out if not using Te, Ti solver
        IF (esrc .NE. 0) THEN
          te(ic)     = f_Te(ic+1)
          ti(ic,ion) = f_Te(ic+1)
        ENDIF

        ne(ic)     = ni(ic,ion)
        pe(ic)     = ne(ic) * te(ic) * ECH
        pi(ic,ion) = ni(ic,ion) * 
     .      (ti(ic,ion) * ECH + mi(ion) * vi(ic,ion)**2)
      ENDDO
 
      RETURN
 99   STOP
      END
c    jrh mod end
c
c ======================================================================
c