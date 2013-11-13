      MODULE SOLF1D

      !constants
      !proton mass in SI units
      DOUBLE PRECISION,PARAMETER::mp=1.6726D-27

      !electron mass in SI units
      DOUBLE PRECISION,PARAMETER::me=9.1094D-31
      !Boltzmann constant in J/eV	
      DOUBLE PRECISION,PARAMETER::B_K=1.6D-19	
      !elementary charge
      DOUBLE PRECISION,PARAMETER::qelement=1.6021892D-19
      !energy transmission factor
      DOUBLE PRECISION::deltae=5.0D0
      !energy transmission factor 
      DOUBLE PRECISION::deltai=3.5D0

      TYPE STRUCT
        DOUBLE PRECISION,DIMENSION(2,2)::A,B,C
        DOUBLE PRECISION,DIMENSION(2)::F		
      END TYPE

      CONTAINS

!------------------------------STEADY-STATE------------------------------
      SUBROUTINE STEADY_STATE(f_n,f_v,f_Te,f_Ti,f_n0,
     .                        f_v0,x_n,x_v,n,Sext_n,
     .                        Sext_v,Sext_Te,Sext_Ti,
     .                        Seir_n,Seir_v,Seir_Te,
     .                        Seir_Ti,Bsol,alphae,
     .                        alphai,beta,NEUT,R,A,IT,
     .                        IT_graph, BC_Dirichlet,
     .                        solve_t)

      IMPLICIT NONE
										
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::f_n,
     .                f_v,f_Te,f_Ti,f_n0,f_v0
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_n,x_v
      INTEGER,INTENT(IN)::n
      !number of grid points (velocity grid)
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Sext_n,
     .             Sext_v,Sext_Te,Sext_Ti,Seir_n,Seir_v,
     .             Seir_Te,Seir_Ti
      !sources
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Bsol
      DOUBLE PRECISION,INTENT(IN)::alphae,alphai,beta
      ! flux limiters	
      INTEGER,INTENT(IN)::NEUT
      !switch for neutrals	
      DOUBLE PRECISION,INTENT(IN)::R
      !recycling coefficient
      INTEGER,INTENT(IN)::A
      !atomic mass number	
      INTEGER,INTENT(IN)::IT,IT_graph,solve_t
      DOUBLE PRECISION,DIMENSION(8),INTENT(IN),OPTIONAL::BC_Dirichlet

      !solution in previous time step
      DOUBLE PRECISION::f0_n(n+1),f0_v(n),f0_Te(n+1),f0_Ti(n+1),
     .                  f0_n0(n+1),f0_v0(n+1)	

      !temperature of neutrals (not an input)
      DOUBLE PRECISION::f_T0(n+1)
	
      !ion mass
      DOUBLE PRECISION::mi,m0
	
      !coefficients of equations
      DOUBLE PRECISION::v_n(n),v_v(n+1),v_Te(n),v_Ti(n),
     .                  v_n0(n),v_v0(n+1) 	
      DOUBLE PRECISION::a_n(n+1),a_n_new(n+1),b_n(n+1),
     .                  Sexp_n(n+1),Simp_n_new(n+1) 
      DOUBLE PRECISION::a_v(n),a_v_new(n),b_v(n),c_v(n),
     .                  c_v_new(n),d_v(n+1),d_v_new(n+1),
     .                  Sexp_v(n),Simp_v_new(n) 
      DOUBLE PRECISION::a_Te(n+1),a_Te_new(n+1),
     .                  b_Te(n+1),c_Te(n+1),c_Te_new(n+1),
     .                  d_Te(n),d_Te_new(n),Sexp_Te(n+1),
     .                  Simp1_Te_new(n+1),Simp2_Te_new(n+1) 
      DOUBLE PRECISION::a_Ti(n+1),a_Ti_new(n+1),b_Ti(n+1),
     .                  c_Ti(n+1),c_Ti_new(n+1),d_Ti(n),
     .                  d_Ti_new(n),Sexp_Ti(n+1),
     .                  Simp1_Ti_new(n+1),Simp2_Ti_new(n+1) 
      DOUBLE PRECISION::a_n0(n+1),a_n0_new(n+1),b_n0(n+1),
     .                  Sexp_n0(n+1),Simp_n0_new(n+1) 
      DOUBLE PRECISION::a_v0(n),a_v0_new(n),b_v0(n),
     .                  c_v0(n),c_v0_new(n),d_v0(n+1),
     .                  d_v0_new(n+1),Sexp_v0(n),
     .                  Simp_v0_new(n) 
      DOUBLE PRECISION::BC_n(2),BC_v(2),BC_T(3,4)
      !boundary conditions arrays
      DOUBLE PRECISION::BC_n0(2),BC_v0(2)
      !boundary conditions arrays	

      !collisions
      DOUBLE PRECISION::ION(n+1),REC(n+1),CX(n+1),EXC(n+1)
      !collision rate coefficients for ionization, 
      !recombination, charge exchange and excitation
	
      !transport
      DOUBLE PRECISION::etai(n+1),kappae(n+1),kappai(n+1),
     .                  lambda(n+1),exchange(n+1)
      !ion viscosity, heat conductivities, coulomb logarithm, 
      !exchange term

      DOUBLE PRECISION::delta_t
      INTEGER::k,i
	
      mi=mp*A
      m0=mi+me
	
      !------------------INITIAL CONDITIONS-------------
      CALL INIT_CONDITIONS(f_n,f_v,f_Te,f_Ti,f_n0,f_v0,
     .                     f_T0,f0_n,f0_v,f0_Te,f0_Ti,
     .                     f0_n0,f0_v0,x_n,x_v,mi,m0,
     .                     B_K,R,N,NEUT)

      !----------------------ITERATION------------------
      DO k=1,IT		

        !----TIMESTEP----				
        delta_t=TIME_STEP(x_n,f_v)				
	
        !---COLLISIONS---		
        CALL ATOM_DATA(ION,REC,CX,EXC,f_n,f_Te,n+1,NEUT)

        !----TRANSPORT---
        CALL TRANSP_COEFF(etai,kappae,kappai,lambda,
     .                    exchange,alphae,alphai,beta,
     .                    x_n,x_v,f_n,f_v,f_Te,f_Ti,Bsol,
     .                    me,mi,B_K,n+1)

        !-----PLASMA-----		
        !coefficients of transport equations		
        CALL COEFFS_N(v_n,a_n,a_n_new,b_n,Sexp_n,
     .                Simp_n_new,x_n,x_v,f_n,f_v,
     .                f_Te,f_n0,n+1,Sext_n,Seir_n,
     .                ION,REC,Bsol)

        CALL COEFFS_V(v_v,a_v,a_v_new,b_v,c_v,
     .                c_v_new,d_v,d_v_new,Sexp_v,
     .                Simp_v_new,x_n,x_v,f_n,f_v,
     .                f_Te,f_Ti,f_n0,f_v0,n,mi,B_K,
     .                Sext_n,Sext_v,Seir_n,Seir_v,ION,
     .                REC,CX,etai,kappai,Bsol)

        CALL COEFFS_T(v_Te,a_Te,a_Te_new,b_Te,c_Te,
     .                c_Te_new,d_Te,d_Te_new,Sexp_Te,
     .                Simp1_Te_new,Simp2_Te_new,v_Ti,
     .                a_Ti ,a_Ti_new,b_Ti,c_Ti,c_Ti_new,
     .                d_Ti,d_Ti_new,Sexp_Ti,
     .                Simp1_Ti_new,Simp2_Ti_new,
     .                x_n,x_v,f_n,f_v,f_Te,f_Ti,f_n0,
     .                f_v0,f_T0,n+1,mi,m0,B_K,Sext_n,
     .                Sext_v,Sext_Te,Sext_Ti,Seir_n,
     .                Seir_v,Seir_Te,Seir_Ti,ION,
     .                REC,CX,EXC,etai,kappae,kappai,
     .                exchange,Bsol)
   
        !boundary conditions		
        CALL BC(BC_n,BC_v,BC_T,x_n,x_v,f_n,f_v,
     .          f_Te,f_Ti,mi,B_K,deltae,deltai,
     .          n+1,etai,kappae,kappai,Bsol,BC_Dirichlet)
  
        !solver
        CALL SOLVER_N(f_n,f0_n,v_n,a_n,a_n_new,b_n,
     .                Sexp_n,Simp_n_new,x_n,x_v,n+1,
     .                delta_t,BC_n)
        CALL SOLVER_V(f_v,f0_v,v_v,a_v,a_v_new,b_v,c_v,
     .                c_v_new,d_v,d_v_new,Sexp_v,
     .                Simp_v_new,x_v,x_n,n,delta_t,BC_v)
        IF (solve_t .EQ. 1) THEN
          CALL SOLVER_T(f_Te,f_Ti,f0_Te,f0_Ti,v_Te,v_Ti,
     .                  a_Te,a_Ti,a_Te_new,a_Ti_new,b_Te,
     .                  b_Ti,c_Te,c_Ti,c_Te_new,c_Ti_new,
     .                  d_Te,d_Ti,d_Te_new,d_Ti_new,
     .                  Sexp_Te,Sexp_Ti,Simp1_Te_new,
     .                  Simp2_Te_new,Simp1_Ti_new,
     .                  Simp2_Ti_new,x_n,x_v,n+1,
     .                  delta_t,BC_T) 
        ENDIF

        !---NEUTRALS----
        IF (NEUT.EQ.1) THEN
          !coefficients of transport equations
          CALL COEFFS_N0(v_n0,a_n0,a_n0_new,b_n0,
     .                   Sexp_n0,Simp_n0_new,f_n,
     .                   f_Te,f_n0,f_v0,n+1,ION,REC)	
          CALL COEFFS_V0(v_v0,a_v0,a_v0_new,b_v0,
     .                   c_v0,c_v0_new,d_v0,d_v0_new,
     .                   Sexp_v0,Simp_v0_new,x_n,x_v,f_n,
     .                   f_v,f_Te,f_n0,f_v0,f_T0,n,mi,m0,
     .                   B_K,REC,CX) 
      
          !boundary conditions
          CALL BC0(BC_n0,BC_v0,f_n,f_v,f_Ti,f_v0,m0,B_K,R,n+1)
      
          !solver  
          CALL SOLVER_N(f_n0,f0_n0,v_n0,a_n0,a_n0_new,
     .                  b_n0,Sexp_n0,Simp_n0_new,x_n,x_v,
     .                  n+1,delta_t,BC_n0)
          CALL SOLVER_V(f_v0,f0_v0,v_v0,a_v0,a_v0_new,
     .                  b_v0,c_v0,c_v0_new,d_v0,d_v0_new,
     .                  Sexp_v0,Simp_v0_new,x_v,x_n,n,delta_t,
     .                  BC_v0)
          f_T0=f_Ti
        ENDIF
      END DO	
	
      !--------------------FINAL OUTPUT-----------------	
c      CALL OUTPUT_FINAL(f_n,f_v,f_Te,f_Ti,f_n0,
c     .                  f_v0,x_n,x_v,Sext_n,Sext_v,Sext_Te,
c     .                  Sext_Ti,Seir_n,Seir_v,Seir_Te,Seir_Ti,
c     .                  Bsol,n)

c      PRINT*,' '
c      PRINT*,'<<<<<<<<<<<<<<<<<<<<<<<<<'
c      PRINT*,'<<    END OF SOLF1D    <<'
c      PRINT*,'<<<<<<<<<<<<<<<<<<<<<<<<<'

      END SUBROUTINE

      !---------------------------INITIAL CONDITIONS---------------------------
      SUBROUTINE INIT_CONDITIONS(f_n,f_v,f_Te,f_Ti,f_n0,
     .                           f_v0,f_T0,f0_n,f0_v,f0_Te,
     .                           f0_Ti,f0_n0,f0_v0,x_n,x_v,
     .                           mi,m0,B_K,R,N,NEUT)
      IMPLICIT NONE
	
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::f_n,f_v,
     .                             f_Te,f_Ti,f_n0,f_v0,f_T0
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::f0_n,f0_v,
     .                                f0_Te,f0_Ti,f0_n0,f0_v0
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_n,x_v
      DOUBLE PRECISION,INTENT(IN)::mi,m0,B_K,R
      !constants	
      INTEGER,INTENT(IN)::N
      !number of velocity points
      INTEGER,INTENT(IN)::NEUT
      !switch for neutral model		
      DOUBLE PRECISION::L
      INTEGER::i

      L=x_n(N+1)-x_n(1)

      !PLASMA
	
      !density
      f_n=1.0D23 
		
      !temperature 
      !f_Te=50.0D0 
      !f_Ti=80.0D0 

      !velocity 	
      f_v=sqrt(B_K*(f_Te(1)+f_Ti(1))/mi)*(x_v)/(0.5D0*L)
      !linear between bc	
    
      !NEUTRALS
	
      !without neutrals
      IF (NEUT.EQ.0) THEN
        f_v0=0.0D0 
        f_n0=0.0D0 
        f_T0=0.0D0
      ENDIF

      !with neutrals	 
      IF (NEUT.EQ.1) THEN       
        f_v0=-sqrt(B_K*f_Ti(1)/m0)*(x_v)/(0.5D0*L) 	        
        f_n0=-R*f_n(1)*f_v(1)/f_v0(1)*((x_n)/(0.5D0*L))**20D0+1.0D16
        f_T0=f_Ti
      ENDIF

      !previous time step
      f0_n=f_n
      f0_v=f_v
      f0_Te=f_Te
      f0_Ti=f_Ti
      f0_n0=f_n0
      f0_v0=f_v0

      END SUBROUTINE

      !--------------------------------SOURCES---------------------------------
      SUBROUTINE SOURCES(Sext_n,Sext_Te,Sext_Ti,x_n,Sn0,
     .                   SEe0,SEi0,sigma,L0)
      IMPLICIT NONE
										
      DOUBLE PRECISION,DIMENSION(:),INTENT(OUT)::Sext_n,
     .                                   Sext_Te,Sext_Ti
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_n
      DOUBLE PRECISION,INTENT(IN)::Sn0,SEe0,SEi0,sigma,L0
	
      Sext_n=Sn0*exp(-(x_n-L0)*(x_n-L0)/(2.0D0*sigma*sigma))
      Sext_Te=SEe0*exp(-(x_n-L0)*(x_n-L0)/(2.0D0*sigma*sigma))
      Sext_Ti=SEi0*exp(-(x_n-L0)*(x_n-L0)/(2.0D0*sigma*sigma))

      END SUBROUTINE

      !--------------------------CONTINUITY EQUATION---------------------------
      !equation: d(Af)/dt + Bd(vf)/dx = Sexp - fSimp 
      !solver: explicit upwind scheme
      SUBROUTINE SOLVER_N(f,f0,v,a,a_new,b,Sexp,Simp_new,
     .                    x_f,x_v,N,delta_t,BC_n)        
      IMPLICIT NONE
										
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::f,f0
      !physical quantity
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::v
      !velocity field	
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::a,a_new,b
      !coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Sexp,Simp_new
      !sources
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_f,x_v
      !grid for function and velocity (staggered grid)
      INTEGER,INTENT(IN)::N
      !number of nodal points of f function
      DOUBLE PRECISION,INTENT(IN)::delta_t
      !time step
      DOUBLE PRECISION,DIMENSION(2),INTENT(IN)::BC_n
      !boundary conditions
      DOUBLE PRECISION::delta_x
      !grid step
      DOUBLE PRECISION::Gamma_up,Gamma_down,Gamma0_up,
     .                  Gamma0_down
      !flux
      DOUBLE PRECISION,DIMENSION(N)::fnew
      !function in next time level
      INTEGER::i

      fnew=f
      DO i=2,N-1
		
        !grid step
        delta_x=x_v(i)-x_v(i-1)
		
        !fluxes
        IF (v(i) .GE. 0.0D0) THEN 
          Gamma_up=f(i)*v(i) 
          Gamma0_up=f0(i)*v(i) 
        ELSE 
          Gamma_up=f(i+1)*v(i) 
          Gamma0_up=f0(i+1)*v(i) 
        ENDIF		
        IF (v(i-1) .GE. 0.0D0) THEN  
          Gamma_down=f(i-1)*v(i-1)
          Gamma0_down=f0(i-1)*v(i-1)
        ELSE
          Gamma_down=f(i)*v(i-1)
          Gamma0_down=f0(i)*v(i-1)
        ENDIF
		
        !next time level
        fnew(i)=(2.0D0*f(i)-0.5D0*f0(i)+delta_t*Sexp(i)
     .           -2.0D0*b(i)*delta_t*(Gamma_up-Gamma_down)
     .           /delta_x+b(i)*delta_t*(Gamma0_up-Gamma0_down)
     .           /delta_x)/(1.5D0+delta_t*Simp_new(i))
	
      END DO
      f0=f
      f=fnew

      IF (BC_n(1).GT.0.0D0) THEN !boundary conditions - dirichlet	
        f(1)=BC_n(1)  
        f(N)=BC_n(2)	
      ELSE                    !boundary conditions - extrapolation
        f(1)=LIN_INTERP(x_f(1),x_f(2),x_f(3),f(2),f(3))
        f(N)=LIN_INTERP(x_f(N),x_f(N-1),x_f(N-2),f(N-1),f(N-2))
      ENDIF	
	
      END SUBROUTINE 

      !---------------------------MOMENTUM EQUATION----------------------------
      !equation: d(Af)/dt + Bd(vf)/dx + Dd(C(df/dx))/dx = Sexp - fSimp 
      !solver: explicit upwind scheme for convective term, C-N for diffusive term 
      SUBROUTINE SOLVER_V(f,f0,v,a,a_new,b,c,c_new,d,d_new,Sexp,Simp_new
     .                    ,x_f,x_v,N,delta_t,BC_v) 

      IMPLICIT NONE
										
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::f,f0
      !physical quantity
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::v
      !velocity field	
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::a,a_new,b,c,c_new,d,
     .                                          d_new
      !coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Sexp,Simp_new
      !sources
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_f,x_v
      !grid for function and velocity (staggered grid)
      INTEGER,INTENT(IN)::N
      !number of nodal points of f function
      DOUBLE PRECISION,INTENT(IN)::delta_t
      !time step
      DOUBLE PRECISION,DIMENSION(2),INTENT(IN)::BC_v
      !boundary conditions
      DOUBLE PRECISION,DIMENSION(N+1,N)::Matrix
      !matrix for algebraic system of equation	
      DOUBLE PRECISION::delta_x,delta_x_1,delta_x_2
      !grid steps
      DOUBLE PRECISION::Gamma_up,Gamma_down,Gamma0_up,Gamma0_down
      !flux
      INTEGER::i	

      !set to zero
      Matrix=0.0D0

      !DISCRETIZATION

      !first point (Dirichlet boundary condition)
      Matrix(1,1)=1.0D0
      Matrix(N+1,1)=BC_v(1) 

      !last point (Dirichlet boundary condition)
      Matrix(N,N)=1.0D0	
      Matrix(N+1,N)=BC_v(2) 

      !middle points
      DO i=2,N-1

        !fluxes
        IF (v(i+1) .GE. 0.0D0) THEN 
          Gamma_up=f(i)*v(i+1) 
          Gamma0_up=f0(i)*v(i+1) 
        ELSE 
          Gamma_up=f(i+1)*v(i+1) 
          Gamma0_up=f0(i+1)*v(i+1) 
        ENDIF		
        IF (v(i) .GE. 0.0D0) THEN
          Gamma_down=f(i-1)*v(i)
          Gamma0_down=f0(i-1)*v(i)
          ELSE
            Gamma_down=f(i)*v(i)
            Gamma0_down=f0(i)*v(i)
          ENDIF

          !grid steps
          delta_x=x_v(i+1)-x_v(i)
          delta_x_1=x_f(i+1)-x_f(i)
          delta_x_2=x_f(i)-x_f(i-1)
		 			
          !coefficients		
          Matrix(i,i)=1.5D0+Simp_new(i)*delta_t-c_new(i)*d_new(i+1)*
     .                delta_t/(delta_x*delta_x_1)-c_new(i)*d_new(i)*
     .                delta_t/(delta_x*delta_x_2)
          Matrix(i+1,i)=delta_t*c_new(i)*d_new(i+1)/
     .                (delta_x*delta_x_1)

          Matrix(i-1,i)=delta_t*c_new(i)*d_new(i)/(delta_x*delta_x_2)

          Matrix(N+1,i)=2.0D0*f(i)-0.5D0*f0(i)+Sexp(i)*delta_t-2.0D0*
     .                  b(i)*delta_t*(Gamma_up-Gamma_down)/delta_x+
     .                  b(i)*delta_t*(Gamma0_up-Gamma0_down)/delta_x

      END DO	

      f0=f

      !SOLVING THE SYSTEM OF EQUATIONS	

      CALL PROGONKA(Matrix,N,f)
	
      END SUBROUTINE

      !-------------------------TEMPERATURE EQUATION---------------------------
      !equation1: d(A1f1)/dt + B1d(v1f1)/dx + D1d(C1(df1/dx))/dx = Sexp1 - f1Simp11 - f2Simp21 
      !equation1: d(A2f2)/dt + B2d(v2f2)/dx + D2d(C2(df2/dx))/dx = Sexp2 - f2Simp12 - f1Simp22
      !solver: explicit upwind scheme for convective term, C-N for diffusive term
      SUBROUTINE SOLVER_T(f1,f2,f01,f02,v1,v2,a1,a2,a1_new,a2_new,b1,
     .                    b2,c1,c2,c1_new,c2_new,d1,d2,d1_new,d2_new,
     .                    Sexp1,Sexp2,Simp11_new,Simp21_new,Simp12_new,
     .                    Simp22_new,x_f,x_v,N,delta_t,BC_T) 
      IMPLICIT NONE
										
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::f1,f2,f01,f02
      !physical quantity
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::v1,v2
      !velocity field	
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::a1,a1_new,b1,c1,c1_new,
     .                                          d1,d1_new
      !coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::a2,a2_new,b2,c2,c2_new,
     .                                          d2,d2_new
      !coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Sexp1,Simp11_new,
     .                                          Simp21_new
      !sources
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Sexp2,Simp12_new,
     .                                          Simp22_new
      !sources
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_f,x_v
      !grid for function and velocity (staggered grid)
      INTEGER,INTENT(IN)::N
      !number of nodal points of f function
      DOUBLE PRECISION,INTENT(IN)::delta_t
      !time step
      DOUBLE PRECISION,DIMENSION(3,4),INTENT(IN)::BC_T
      !boundary conditions
      DOUBLE PRECISION::delta_x,delta_x_1,delta_x_2
      !grid steps
      DOUBLE PRECISION::Gamma1_up,Gamma1_down,Gamma2_up,Gamma2_down
      !flux
      DOUBLE PRECISION::Gamma01_up,Gamma01_down,Gamma02_up,Gamma02_down
      !flux
      TYPE(STRUCT),DIMENSION(N)::coeffs
      !structure of coefficients for MATRIX_PROGONKA
      INTEGER::i		
	

      DO i=1,N 
        coeffs(i)%A=0D0
        coeffs(i)%A=0D0
        coeffs(i)%B=0D0
        coeffs(i)%C=0D0
        coeffs(i)%F=0D0
      END DO
	
      !DISCRETIZATION

      !boundary conditions
      coeffs(1)%B(1,1)=BC_T(2,1)-BC_T(1,1)/(x_f(2)-x_f(1))
      !BC_T(1,1)*df1/dx + BC_T(2,1)*f1 = BC_T(3,1) (BC for f1 at point 1)
      coeffs(1)%B(2,2)=BC_T(2,2)-BC_T(1,2)/(x_f(2)-x_f(1))
      !BC_T(1,2)*df2/dx + BC_T(2,2)*f2 = BC_T(3,2) (BC for f2 at point 1)
      coeffs(1)%C(1,1)=BC_T(1,1)/(x_f(2)-x_f(1))
      !BC_T(1,3)*df1/dx + BC_T(2,3)*f1 = BC_T(3,3) (BC for f1 at point N)
      coeffs(1)%C(2,2)=BC_T(1,2)/(x_f(2)-x_f(1))
      !BC_T(1,4)*df2/dx + BC_T(2,4)*f2 = BC_T(3,4) (BC for f2 at point N)
      coeffs(1)%F(1)=BC_T(3,1)
      coeffs(1)%F(2)=BC_T(3,2)
      coeffs(N)%A(1,1)=-BC_T(1,3)/(x_f(N)-x_f(N-1))	
      coeffs(N)%A(2,2)=-BC_T(1,4)/(x_f(N)-x_f(N-1))
      coeffs(N)%B(1,1)=BC_T(2,3)+BC_T(1,3)/(x_f(N)-x_f(N-1))
      coeffs(N)%B(2,2)=BC_T(2,4)+BC_T(1,4)/(x_f(N)-x_f(N-1))
      coeffs(N)%F(1)=BC_T(3,3)
      coeffs(N)%F(2)=BC_T(3,4)

      !middle points 	
      DO i=2,N-1
		
        !fluxes
        IF (v1(i) .GE. 0.0D0) THEN 
          Gamma1_up=f1(i)*v1(i) 
          Gamma01_up=f01(i)*v1(i) 
        ELSE 
          Gamma1_up=f1(i+1)*v1(i) 
          Gamma01_up=f01(i+1)*v1(i) 
        ENDIF		
        IF (v1(i-1) .GE. 0.0D0) THEN
          Gamma1_down=f1(i-1)*v1(i-1)
          Gamma01_down=f01(i-1)*v1(i-1)
        ELSE
          Gamma1_down=f1(i)*v1(i-1)
          Gamma01_down=f01(i)*v1(i-1)
        ENDIF
        IF (v2(i) .GE. 0.0D0) THEN 
          Gamma2_up=f2(i)*v2(i) 
          Gamma02_up=f02(i)*v2(i) 
        ELSE 
          Gamma2_up=f2(i+1)*v2(i) 
          Gamma02_up=f02(i+1)*v2(i) 
        ENDIF		
        IF (v2(i-1) .GE. 0.0D0) THEN
          Gamma2_down=f2(i-1)*v2(i-1)
          Gamma02_down=f02(i-1)*v2(i-1)
        ELSE
          Gamma2_down=f2(i)*v2(i-1)
          Gamma02_down=f02(i)*v2(i-1)
        ENDIF
	
		!grid step
        delta_x=x_v(i)-x_v(i-1)
        delta_x_1=x_f(i+1)-x_f(i)
        delta_x_2=x_f(i)-x_f(i-1)

        !coefficients		
        coeffs(i)%B(1,1)=1.5D0-delta_t*c1_new(i)*d1_new(i)/
     .                   (delta_x*delta_x_1)-delta_t*c1_new(i)*
     .                   d1_new(i-1)/(delta_x*delta_x_2)+
     .                   delta_t*Simp11_new(i)
        coeffs(i)%A(1,1)=delta_t*c1_new(i)*d1_new(i-1)/
     .                   (delta_x*delta_x_2)
        coeffs(i)%C(1,1)=delta_t*c1_new(i)*d1_new(i)/(delta_x*delta_x_1)
        coeffs(i)%F(1)=2.0D0*f1(i)-0.5D0*f01(i)+Sexp1(i)*
     .                 delta_t-2.0D0*b1(i)*delta_t*
     .                 (Gamma1_up-Gamma1_down)/delta_x+b1(i)*delta_t*
     .                 (Gamma01_up-Gamma01_down)/delta_x
        coeffs(i)%B(2,1)=delta_t*Simp21_new(i)

        coeffs(i)%B(2,2)=1.5D0-delta_t*c2_new(i)*d2_new(i)/
     .                   (delta_x*delta_x_1)-delta_t*c2_new(i)*
     .                   d2_new(i-1)/(delta_x*delta_x_2)+delta_t*
     .                   Simp12_new(i)
        coeffs(i)%A(2,2)=delta_t*c2_new(i)*d2_new(i-1)/
     .                   (delta_x*delta_x_2)		
        coeffs(i)%C(2,2)=delta_t*c2_new(i)*d2_new(i)/(delta_x*delta_x_1)
        coeffs(i)%F(2)=2.0D0*f2(i)-0.5D0*f02(i)+Sexp2(i)*delta_t-2.0D0*
     .                 b2(i)*delta_t*(Gamma2_up-Gamma2_down)/
     .                 delta_x+b2(i)*delta_t*(Gamma02_up-Gamma02_down)/
     .                 delta_x
        coeffs(i)%B(1,2)=delta_t*Simp22_new(i)

      END DO

      f01=f1
      f02=f2

      !SOLVING THE SYSTEM OF EQUATIONS
		
      CALL MATRIX_PROGONKA(coeffs,N,f1,f2)

      END SUBROUTINE

      !-------------------------BOUNDARY CONDITIONS----------------------------
      SUBROUTINE BC(BC_n,BC_v,BC_T,x_n,x_v,f_n,f_v,f_Te,f_Ti,mi,B_K,
     .              deltae,deltai,N,etai,kappae,kappai,Bsol,
     .              BC_Dirichlet)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(2),INTENT(INOUT)::BC_n,BC_v
      !boundary conditions for density and velocity
      DOUBLE PRECISION,DIMENSION(3,4),INTENT(INOUT)::BC_T
      !boundary conditions for temperature
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_n,x_v
      !grid
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::f_n,f_v,f_Te,f_Ti
      !density, velocity and temperature
      DOUBLE PRECISION,INTENT(IN)::mi,B_K
      !constants
      DOUBLE PRECISION,INTENT(IN)::deltae,deltai
      !heat transmission coefficients
      INTEGER,INTENT(IN)::N
      !number of temperature points	
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::etai,kappae,kappai
      !transport coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Bsol
      !magnetic field
      DOUBLE PRECISION::sBsol(N-1)
      !staggered functions		
      DOUBLE PRECISION,DIMENSION(8),INTENT(IN),OPTIONAL::BC_Dirichlet
      DOUBLE PRECISION::cs_1,cs_N,kappae_1,kappae_N,kappai_1,kappai_N
      DOUBLE PRECISION::etai_1,etai_N,etaterm_1,etaterm_N
	
      !staggered functions
      sBsol=STAGGER1(Bsol,N,x_n,x_v)

      !quantities at the boundaries
      cs_1=-sqrt(B_K*(f_Te(1)+f_Ti(1))/mi) 	!sound speed
      cs_N=sqrt(B_K*(f_Te(N)+f_Ti(N))/mi)	!sound speed

      kappae_1=kappae(1) 			!thermal conductivity
      kappae_N=kappae(N)			!thermal conductivity
      kappai_1=kappai(1) 			!thermal conductivity
      kappai_N=kappai(N)  			!thermal conductivity
      etai_1=etai(1)				!ion viscosity
      etai_N=etai(N)				!ion viscosity
	
      etaterm_1=-etai_1*((f_v(2)-f_v(1))/(x_v(2)-x_v(1)))/(B_K*f_n(1))
      etaterm_N=-etai_N*((f_v(N-1)-f_v(N-2))/(x_v(N-1)-x_v(N-2)))/
     .          (B_K*f_n(N))

      !PLASMA DENSITY
      IF (PRESENT(BC_Dirichlet)) THEN 
        !Dirichlet
        BC_n(1)=BC_Dirichlet(1)
        BC_n(2)=BC_Dirichlet(2)
      ELSE	
        BC_n(1)=-1
        BC_n(2)=-1 
      ENDIF


      !PLASMA VELOCITY
      IF (PRESENT(BC_Dirichlet)) THEN 
        !Dirichlet
        BC_v(1)=BC_Dirichlet(7)
        BC_v(2)=BC_Dirichlet(8)
      ELSE	
        BC_v(1)=cs_1 
        BC_v(2)=cs_N 
      ENDIF
	
      !PLASMA TEMPERATURE
      !BC_T(1,1)*df1/dx + BC_T(2,1)*f1 = BC_T(3,1) (BC for f1 at point 1)
      !BC_T(1,2)*df2/dx + BC_T(2,2)*f2 = BC_T(3,2) (BC for f2 at point 1)
      !BC_T(1,3)*df1/dx + BC_T(2,3)*f1 = BC_T(3,3) (BC for f1 at point N)
      !BC_T(1,4)*df2/dx + BC_T(2,4)*f2 = BC_T(3,4) (BC for f2 at point N)

      IF (PRESENT(BC_Dirichlet)) THEN 
        !Dirichlet (f1=1, f2=1)
        BC_T(1,1)=0; BC_T(2,1)=1; BC_T(3,1)=BC_Dirichlet(3)
        BC_T(1,2)=0; BC_T(2,2)=1; BC_T(3,2)=BC_Dirichlet(5)
        BC_T(1,3)=0; BC_T(2,3)=1; BC_T(3,3)=BC_Dirichlet(4)
        BC_T(1,4)=0; BC_T(2,4)=1; BC_T(3,4)=BC_Dirichlet(6)
      ELSE	
        !Newton
        BC_T(1,1)=kappae_1/(f_n(1)*cs_1)
        BC_T(2,1)=(deltae-2.5D0)
        BC_T(3,1)=0D0  
        BC_T(1,2)=kappai_1/(cs_1*f_n(1))
        BC_T(2,2)=(deltai-2.5D0-0.5D0)
        BC_T(3,2)=0.5D0*f_Te(1)+etaterm_1
        BC_T(1,3)=kappae_N/(f_n(N)*cs_N)
        BC_T(2,3)=(deltae-2.5D0)
        BC_T(3,3)=0D0
        BC_T(1,4)=kappai_N/(cs_N*f_n(N))
        BC_T(2,4)=(deltai-2.5D0-0.5D0)
        BC_T(3,4)=0.5D0*f_Te(N)+etaterm_N
      ENDIF

      !Neumann (df1/dx=0, df2/dx=0)
      !BC_T(1,1)=1
      !BC_T(2,1)=0
      !BC_T(3,1)=0
      !BC_T(1,2)=1
      !BC_T(2,2)=0
      !BC_T(3,2)=0
      !BC_T(1,3)=1
      !BC_T(2,3)=0
      !BC_T(3,3)=0
      !BC_T(1,4)=1
      !BC_T(2,4)=0
      !BC_T(3,4)=0

      !Newton
      !BC_T(1,2)=kappai_1/(cs_1*f_n(1))
      !BC_T(2,2)=(deltai-2.5-0.5-0.5*f_Te(1)/f_Ti(1)-etaterm_1/f_Ti(1))
      !BC_T(3,2)=0	
      !BC_T(1,4)=kappai_N/(cs_N*f_n(N))
      !BC_T(2,4)=(deltai-2.5-0.5-0.5*f_Te(N)/f_Ti(N)-etaterm_N/f_Ti(N))
      !BC_T(3,4)=0
      !IF (BC_T(2,2).LE.0.0) BC_T(2,2)=0.0 !gradient limit
      !IF (BC_T(2,4).LE.0.0) BC_T(2,4)=0.0
	
      !Newton equivalent to SOLPS5.0
      !BC_T(1,1)=kappae_1/(f_n(1)*cs_1)
      !BC_T(2,1)=(deltae-2.5D0)
      !BC_T(3,1)=0D0  
      !BC_T(1,2)=kappai_1/(cs_1*f_n(1))
      !BC_T(2,2)=(deltai-2.5D0)
      !BC_T(3,2)=0D0
      !BC_T(1,3)=kappae_N/(f_n(N)*cs_N)
      !BC_T(2,3)=(deltae-2.5D0)
      !BC_T(3,3)=0D0
      !BC_T(1,4)=kappai_N/(cs_N*f_n(N))
      !BC_T(2,4)=(deltai-2.5D0)
      !BC_T(3,4)=0D0
		
      END SUBROUTINE

      !---------------------NEUTRAL BOUNDARY CONDITIONS------------------------
      SUBROUTINE BC0(BC_n0,BC_v0,f_n,f_v,f_Ti,f_v0,m0,B_K,R,N)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(2),INTENT(INOUT)::BC_n0,BC_v0
      !boundary conditions for density and velocity    
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::f_n,f_v,f_Ti,f_v0
      !plasma density and velocity,ion temperature, neutral velocity
      DOUBLE PRECISION,INTENT(IN)::m0,B_K,R
      !constants and recycling coefficient	
      INTEGER,INTENT(IN)::N
      !number of density points
	
      !NEUTRAL VELOCITY
      BC_v0(1)=sqrt(B_K*f_Ti(1)/m0)
      BC_v0(2)=-sqrt(B_K*f_Ti(N)/m0) 

      !NEUTRAL DENSITY
      BC_n0(1)=-R*f_n(1)*f_v(1)/f_v0(1)
      BC_n0(2)=-R*f_n(N)*f_v(N-1)/f_v0(N-1)			
		
      END SUBROUTINE

      !-------------------------------PROGONKA---------------------------------
      SUBROUTINE PROGONKA(Matrix,N,x) 
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN)::Matrix
      !Matrix*x=b, b is last column of Matrix
      INTEGER,INTENT(IN)::N
      !size of matrix N+1 x N
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::x
      !solution 
      DOUBLE PRECISION,DIMENSION(N)::alpha,beta
      !coefficients
      INTEGER::i	

      !boundary condition
      alpha(2)=0D0
      beta(2)=Matrix(N+1,1)

      !coefficients alpha and beta
      DO i=2,N-1
        alpha(i+1)=-Matrix(i+1,i)/(Matrix(i-1,i)*alpha(i)+Matrix(i,i))
        beta(i+1)=(Matrix(N+1,i)-Matrix(i-1,i)*beta(i))/(Matrix(i-1,i)*
     .            alpha(i)+Matrix(i,i))
      END DO

      !boundary condition
      x(N)=Matrix(N+1,N)

      !solution
      DO i=N-1,1,-1
        x(i)=alpha(i+1)*x(i+1)+beta(i+1)
      END DO

      END SUBROUTINE

      !---------------------------MATRIX PROGONKA------------------------------
      SUBROUTINE MATRIX_PROGONKA(coeffs,N,f1,f2) 
      IMPLICIT NONE
      INTEGER,INTENT(IN)::N                               !size of solution
      TYPE(STRUCT),DIMENSION(:),INTENT(IN)::coeffs        !input coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::f1,f2  !solution
      DOUBLE PRECISION,DIMENSION(2,2)::INVERTED           
      DOUBLE PRECISION,DIMENSION(2)::f
      INTEGER::i

      TYPE STRUCT1
        DOUBLE PRECISION,DIMENSION(2,2)::ALPHA
        DOUBLE PRECISION,DIMENSION(2)::BETA		
      END TYPE

      TYPE(STRUCT1),DIMENSION(N)::MP

      !ALPHA2 and BETA2 from boundary condition
      INVERTED=INVERT(coeffs(1)%B)
      MP(2)%ALPHA=-MULTIPLY(INVERTED,coeffs(1)%C)
      MP(2)%BETA=MULTIPLY_VEC(INVERTED,coeffs(1)%F)

      !ALPHA3 and BETA3 to ALPHAN and BETAN
      DO i=2,N-1
        INVERTED=INVERT(MULTIPLY(coeffs(i)%A,MP(i)%ALPHA)+coeffs(i)%B)
        MP(i+1)%ALPHA=-MULTIPLY(INVERTED,coeffs(i)%C)
        MP(i+1)%BETA=MULTIPLY_VEC(INVERTED,(coeffs(i)%F-
     .               MULTIPLY_VEC(coeffs(i)%A,MP(i)%BETA)))
      END DO

      !fN from boundary condition
      INVERTED=INVERT(MULTIPLY(coeffs(N)%A,MP(N)%ALPHA)+coeffs(N)%B)
      f=MULTIPLY_VEC(INVERTED,coeffs(N)%F-
     .               MULTIPLY_VEC(coeffs(N)%A,MP(N)%BETA)) 
      f1(N)=f(1)
      f2(N)=f(2)

      !fN-1 to f1
      DO i=N-1,1,-1
        f=MULTIPLY_VEC(MP(i+1)%ALPHA,f)+MP(i+1)%BETA
        f1(i)=f(1)
        f2(i)=f(2)
      END DO	
	
      END SUBROUTINE

      !-------------------------------INVERT-----------------------------------
      FUNCTION INVERT(M) !inverse of 2x2 array
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(2,2)::INVERT
      DOUBLE PRECISION,DIMENSION(2,2),INTENT(IN)::M
      DOUBLE PRECISION::K,L

      K=M(1,2)/(M(2,1)*M(1,2)-M(1,1)*M(2,2))
      L=M(1,1)/(-M(2,1)*M(1,2)+M(1,1)*M(2,2))
      INVERT(2,2)=L
      INVERT(1,2)=K
      INVERT(1,1)=(1.0D0-K*M(2,1))/M(1,1)
      INVERT(2,1)=-L*M(2,1)/M(1,1)
      END FUNCTION

      !-----------------------------MULTIPLY-----------------------------------
      FUNCTION MULTIPLY(M,N) !multiplication of 2x2 arrays
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(2,2)::MULTIPLY
      DOUBLE PRECISION,DIMENSION(2,2),INTENT(IN)::M,N
	
      MULTIPLY(1,1)=M(1,1)*N(1,1)+M(2,1)*N(1,2)
      MULTIPLY(2,1)=M(1,1)*N(2,1)+M(2,1)*N(2,2)
      MULTIPLY(1,2)=M(1,2)*N(1,1)+M(2,2)*N(1,2)
      MULTIPLY(2,2)=M(1,2)*N(2,1)+M(2,2)*N(2,2)
      END FUNCTION

      !-----------------------------MULTIPLY-----------------------------------
      FUNCTION MULTIPLY_VEC(M,vec) !multiplication of vector and 2x2 array
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(2)::MULTIPLY_VEC
      DOUBLE PRECISION,DIMENSION(2,2),INTENT(IN)::M
      DOUBLE PRECISION,DIMENSION(2),INTENT(IN)::vec
	
      MULTIPLY_VEC(1)=M(1,1)*VEC(1)+M(2,1)*VEC(2)    
      MULTIPLY_VEC(2)=M(1,2)*VEC(1)+M(2,2)*VEC(2) 
      END FUNCTION

      !-------------------------COULOMB LOGARITHM------------------------------
      SUBROUTINE TRANSP_COEFF(etai,kappae,kappai,lambda,exchange,alphae,
     .                        alphai,beta,x_n,x_v,f_n,f_v,f_Te,f_Ti,
     .                        Bsol,me,mi,B_K,N)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::etai,kappae,kappai,
     .                              lambda,exchange
      !transport coefficients
      DOUBLE PRECISION,INTENT(IN)::alphae,alphai,beta
      !flux limiters
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_n,x_v
      !grid
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::f_n,f_v,f_Te,f_Ti
      !density, velocity, temperature
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Bsol
      !magnetic field
      DOUBLE PRECISION,INTENT(IN)::me,mi,B_K
      !constants
      INTEGER,INTENT(IN)::N
      !number of points of density grid
      DOUBLE PRECISION::der_Te(N-1),der_Ti(N-1),sder_Te(N),sder_Ti(N),
     .                  der_v(N)
      !gradients
      DOUBLE PRECISION::s_v(N)
      !staggered functions
      INTEGER::i

      !---TRANSPORT COEFFICIENTS
	
      !coulomb logarithm
      lambda=17.0D0
	!DO i=1,N
		!NRL formulary
		!IF (f_Te(i).LT.(f_Ti(i)*me/mi)) THEN 
		!	lambda(i)=30.0D0-LOG(SQRT(f_n(i)*1.0D-6)*(f_Ti(i)**(-1.5D0)))
		!ENDIF
	
		!IF ((f_Te(i).GE.(f_Ti(i)*me/mi)).AND.(f_Te(i).GT.10.0D0)) THEN
		!	lambda(i)=24.0D0-LOG(SQRT(f_n(i)*1.0D-6)/f_Te(i))
		!ENDIF
	
		!IF ((f_Te(i).GE.(f_Ti(i)*me/mi)).AND.(f_Te(i).LE.10.0D0)) THEN
		!	lambda(i)=23.0D0-LOG(SQRT(f_n(i)*1.0D-6)*(f_Te(i)**(-1.5D0)))
		!ENDIF	

		!Wesson
		!lambda(i)=15.2D0-0.5D0*LOG(f_n(i)*1.0D-20)+LOG(f_Te(i)*1.0D-3)			
	!END DO	
	
      !ion viscosity
      etai=4.9D26*B_K*sqrt(mi)*(f_Ti**(2.5D0))/lambda 

      !electron heat conductivity
      kappae=1.1D12*B_K*(f_Te**(2.5D0))/(lambda*me)
	
      !ion heat conductivity
      kappai=19.9D26*B_K*sqrt(mi)*(f_Ti**(2.5D0))/(lambda*mi)

      !exchange term
      exchange=2.0D0*me*f_n*lambda/(0.344D12*mi*(f_Te**(1.5D0))) 

      !---FLUX LIMITERS
		
      !gradients for flux limited coefficients
      DO i=2,N-1 
        der_v(i)=(f_v(i)-f_v(i-1))/(x_v(i)-x_v(i-1))	
      END DO
      der_v(1)=(f_v(2)-f_v(1))/(x_v(2)-x_v(1))
      der_v(N)=(f_v(N-1)-f_v(N-2))/(x_v(N-1)-x_v(N-2))

      DO i=1,N-1 
        der_Te(i)=(f_Te(i+1)-f_Te(i))/(x_n(i+1)-x_n(i))
        der_Ti(i)=(f_Ti(i+1)-f_Ti(i))/(x_n(i+1)-x_n(i))		
      END DO
      sder_Te=STAGGER2(der_Te,N-1,x_n,x_v)
      sder_Ti=STAGGER2(der_Ti,N-1,x_n,x_v)
		
      s_v=STAGGER2(f_v,N-1,x_n,x_v)

      !flux limited coefficients
      etai=etai/(1.0D0+etai*abs(der_v)/(beta*f_n*B_K*f_Ti))		
      kappae=kappae/(1.0D0+kappae*abs(sder_Te)/
     .       (alphae*f_n*f_Te*sqrt(B_K*f_Te/me)))
      kappai=kappai/(1.0D0+kappai*abs(sder_Ti)/
     .       (alphai*f_n*f_Ti*sqrt(B_K*f_Ti/mi)))

      END SUBROUTINE

      !------------------------------HEAT FLUX---------------------------------
      SUBROUTINE HEAT_FLUX(heatfluxe,heatfluxi,thermalhfe,thermalhfi,
     .                     f_n,f_v,f_Te,f_Ti,x_n,x_v,N,mi,B_K,etai,
     .                     kappae,kappai,Bsol)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::heatfluxe,heatfluxi,
     .                              thermalhfe,thermalhfi
      !total heat flux and thermal heat flux
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::f_n,f_v,f_Te,f_Ti
      !density, velocity and temperature
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_n,x_v
      !grid
      INTEGER,INTENT(IN)::N
      !number of density points
      DOUBLE PRECISION,INTENT(IN)::mi,B_K
      !constants	
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::etai,kappae,kappai
      !transport coefficients				 
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Bsol
      !magnetic field
      DOUBLE PRECISION,DIMENSION(N)::qe,qi
      !total heat flux
      DOUBLE PRECISION::s_Te(N-1),s_Ti(N-1),s_v(N),sBsol(N-1)
      !staggered functions		
      INTEGER::i

      !staggered functions
      s_Te=STAGGER1(f_Te,N,x_n,x_v)
      s_Ti=STAGGER1(f_Ti,N,x_n,x_v)
      s_v=STAGGER2(f_v,N-1,x_n,x_v)
      sBsol=STAGGER1(Bsol,N,x_n,x_v)
	
      !kinetic terms
      qe=2.5D0*f_n*B_K*f_Te*s_v
      qi=2.5D0*f_n*B_K*f_Ti*s_v+0.5D0*mi*f_n*(s_v)**3

      !thermal flux and viscosity term
      DO i=2,N-1
        qe(i)=qe(i)-kappae(i)*B_K*(s_Te(i)-s_Te(i-1))/(x_v(i)-x_v(i-1))
        qi(i)=qi(i)-kappai(i)*B_K*(s_Ti(i)-s_Ti(i-1))/(x_v(i)-x_v(i-1))
        qi(i)=qi(i)-s_v(i)*etai(i)*1.0D0/SQRT(Bsol(i))*(sqrt(sBsol(i))*
     .        f_v(i)-sqrt(sBsol(i-1))*f_v(i-1))/(x_v(i)-x_v(i-1))
        thermalhfe(i)=-kappae(i)*B_K*(s_Te(i)-s_Te(i-1))/
     .                (x_v(i)-x_v(i-1))
        thermalhfi(i)=-kappai(i)*B_K*(s_Ti(i)-s_Ti(i-1))/
     .                (x_v(i)-x_v(i-1))
      END DO	
      qe(1)=qe(1)-kappae(1)*B_K*(f_Te(2)-f_Te(1))/(x_n(2)-x_n(1))
      qi(1)=qi(1)-kappai(1)*B_K*(f_Ti(2)-f_Ti(1))/(x_n(2)-x_n(1))
	
      qi(1)=qi(1)-s_v(1)*etai(1)*1.0D0/SQRT(Bsol(1))*
     .      (sqrt(Bsol(2))*s_v(2)-sqrt(Bsol(1))*s_v(1))/(x_n(2)-x_n(1))
	
      qe(N)=qe(N)-kappae(N)*B_K*(f_Te(N)-f_Te(N-1))/(x_n(N)-x_n(N-1))
      qi(N)=qi(N)-kappai(N)*B_K*(f_Ti(N)-f_Ti(N-1))/(x_n(N)-x_n(N-1))
	
      qi(N)=qi(N)-s_v(N)*etai(N)*1.0D0/SQRT(Bsol(N))*
     .      (sqrt(Bsol(N))*s_v(N)-sqrt(Bsol(N-1))*s_v(N-1))/
     .      (x_n(N)-x_n(N-1))
	
      thermalhfe(1)=-kappae(1)*B_K*(f_Te(2)-f_Te(1))/(x_n(2)-x_n(1))
      thermalhfi(1)=-kappai(1)*B_K*(f_Ti(2)-f_Ti(1))/(x_n(2)-x_n(1))
      thermalhfe(N)=-kappae(N)*B_K*(f_Te(N)-f_Te(N-1))/(x_n(N)-x_n(N-1))
      thermalhfi(N)=-kappai(N)*B_K*(f_Ti(N)-f_Ti(N-1))/(x_n(N)-x_n(N-1))

      heatfluxe=qe
      heatfluxi=qi

      END SUBROUTINE

      !----------------------------ELECTRIC FIELD------------------------------
      SUBROUTINE E_FIELD(Efield,f_n,f_Te,x_n,x_v,N,B_K,qelement)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::Efield	!electric field
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::f_n,f_Te	!density and temperature
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_n,x_v	!grid
      INTEGER,INTENT(IN)::N					!number of density points
      DOUBLE PRECISION,INTENT(IN)::B_K,qelement		!constants	
      DOUBLE PRECISION::s_n(N-1),s_Te(N-1)			!staggered functions
      DOUBLE PRECISION::derivation1,derivation2
      INTEGER::i

      !staggered functions
      s_n=STAGGER1(f_n,N,x_n,x_v)
      s_Te=STAGGER1(f_Te,N,x_n,x_v)

      !electric field
      DO i=2,N-1
        derivation1=(s_Te(i)-s_Te(i-1))/(x_v(i)-x_v(i-1))
        derivation2=(s_n(i)*s_Te(i)-s_n(i-1)*s_Te(i-1))/
     .              (x_v(i)-x_v(i-1))
        Efield(i)=-B_K*(0.71*derivation1+derivation2/f_n(i))/qelement
      END DO
      derivation1=(f_Te(2)-f_Te(1))/(x_n(2)-x_n(1))
      derivation2=(f_Te(2)*f_n(2)-f_Te(1)*f_n(1))/(x_n(2)-x_n(1))	
      Efield(1)=-B_K*(0.71*derivation1+derivation2/f_n(1))/qelement
      derivation1=(f_Te(N)-f_Te(N-1))/(x_n(N)-x_n(N-1))
      derivation2=(f_Te(N)*f_n(N)-f_Te(N-1)*f_n(N-1))/(x_n(N)-x_n(N-1))
      Efield(N)=-B_K*(0.71*derivation1+derivation2/f_n(N))/qelement

      END SUBROUTINE
	
      !-----------------COEFFICIENTS FOR CONTINUITY EQUATION-------------------
      SUBROUTINE COEFFS_N(v_n,a_n,a_n_new,b_n,Sexp_n,Simp_n_new,x_n,x_v,
     .                    f_n,f_v,f_Te,f_n0,N,Sext_n,Seir_n,ION,REC,
     .                    Bsol)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::v_n
      !velocity field	
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::a_n,a_n_new,b_n
      !coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::Sexp_n,Simp_n_new
      !sources
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_n,x_v
      !grid
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::f_n,f_v,f_Te,f_n0
      !density, velocity, electron temperature and neutral density
      INTEGER,INTENT(IN)::N
      !number of density points
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Sext_n,Seir_n
      !external source of particles
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::ION,REC
      !collision rate coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Bsol
      !magnetic field
      DOUBLE PRECISION::sBsol(N-1)
      !staggered functions

      !staggered functions
      sBsol=STAGGER1(Bsol,N,x_n,x_v)

      !velocity field
      v_n=f_v/sBsol

      !coefficients
      a_n=1.0D0
      a_n_new=a_n
      b_n=Bsol

      !sources	
      Sexp_n=Sext_n+Seir_n
      Simp_n_new=f_n*REC-f_n0*ION
		
      END SUBROUTINE

      !------------------COEFFICIENTS FOR MOMENTUM EQUATION--------------------
      SUBROUTINE COEFFS_V(v_v,a_v,a_v_new,b_v,c_v,c_v_new,d_v,d_v_new,
     .                    Sexp_v,Simp_v_new,x_n,x_v,f_n,f_v,f_Te,f_Ti,
     .                    f_n0,f_v0,N,mi,B_K,Sext_n,Sext_v,Seir_n,
     .                    Seir_v,ION,REC,CX,etai,kappai,Bsol)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::v_v
      !velocity field	
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::a_v,a_v_new,b_v,
     .                                             c_v,c_v_new,d_v,
     .                                             d_v_new
      !coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::Sexp_v,Simp_v_new
      !sources
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_n,x_v
      !grid
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::f_n,f_v,f_Te,f_Ti,f_n0,
     .                                          f_v0
      !density, velocity and temperature and density and velocity of neutrals	
      INTEGER,INTENT(IN)::N
      !number of velocity points
      DOUBLE PRECISION,INTENT(IN)::mi,B_K
      !constants	
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Sext_n,Sext_v,Seir_n,
     .                                          Seir_v
      !sources of particles and momentum
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::ION,REC,CX
      !collision rate coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::etai,kappai
      !transport coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Bsol
      !magnetic field
      DOUBLE PRECISION::sION(N),sREC(N),sCX(N),s_v(N+1),s_n(N),s_n0(N),
     .                  sSext(N),sSeirn(N),sSeirv(N),sBsol(N),setai(N),
     .                  sSextv(N)
      !staggered functions
      DOUBLE PRECISION::dBdx_v(N),dBdx_n(N+1),dpdx_v(N),dndx_v(N)
      !derivatives dB/dx, dp/dx, dn/dx
      INTEGER::i

      !B grad and derivatives
      DO i=1,N
        dBdx_v(i)=(Bsol(i+1)-Bsol(i))/(x_n(i+1)-x_n(i))
        dpdx_v(i)=B_K*(f_n(i+1)*(f_Te(i+1)+f_Ti(i+1))-f_n(i)*
     .            (f_Te(i)+f_Ti(i)))/(x_n(i+1)-x_n(i))
        dndx_v(i)=(f_n(i+1)-f_n(i))/(x_n(i+1)-x_n(i))		
      END DO
	
      !staggered functions
      sION=STAGGER1(ION,N+1,x_n,x_v)
      sREC=STAGGER1(REC,N+1,x_n,x_v)
      sCX=STAGGER1(CX,N+1,x_n,x_v)  
      sSext=STAGGER1(Sext_n,N+1,x_n,x_v)
      sSeirn=STAGGER1(Seir_n,N+1,x_n,x_v)
      sSeirv=STAGGER1(Seir_v,N+1,x_n,x_v)
      sSextv=STAGGER1(Sext_v,N+1,x_n,x_v)
      s_n=STAGGER1(f_n,N+1,x_n,x_v)
      s_n0=STAGGER1(f_n0,N+1,x_n,x_v)
      s_v=STAGGER2(f_v,N,x_n,x_v)
      sBsol=STAGGER1(Bsol,N+1,x_n,x_v)
      dBdx_n=STAGGER2(dBdx_v,N,x_n,x_v)
      setai=STAGGER1(etai,N+1,x_n,x_v)

      !velocity field (in midpoints of velocity grid)
      v_v=0.5D0*s_v-etai*dBdx_n/(2.0D0*Bsol*mi*f_n)

      !coefficients
      a_v=1.0D0
      a_v_new=a_v
      b_v=1.0D0
      c_v=-sBsol**(1.5)/(mi*s_n)
      c_v_new=c_v
      d_v=etai*Bsol**(-1.5)
      d_v_new=d_v

      !sources
      Sexp_v=-dpdx_v/(mi*s_n)
      Sexp_v=Sexp_v+f_v*setai*dBdx_v*dndx_v/(2.0D0*sBsol*mi*s_n*s_n)
      Sexp_v=Sexp_v-3.0D0*f_v*setai*dBdx_v*dBdx_v/
     .       (4.0D0*mi*s_n*sBsol*sBsol)
      Sexp_v=Sexp_v+(sSeirv+sSextv)/(mi*s_n)+s_n0*f_v0*(sION+sCX)
      Simp_v_new=(sSext+sSeirn)/s_n+s_n0*(sION+sCX)

      END SUBROUTINE

      !----------------COEFFICIENTS FOR TEMPERATURE EQUATIONS------------------
      !equation1: d(A1f1)/dt + B1d(v1f1)/dx + D1d(C1(df1/dx))/dx = Sexp1 - f1Simp11 - f2Simp21 
      !equation2: d(A2f2)/dt + B2d(v2f2)/dx + D2d(C2(df2/dx))/dx = Sexp2 - f2Simp12 - f1Simp22
      SUBROUTINE COEFFS_T(v_Te,a_Te,a_Te_new,b_Te,c_Te,c_Te_new,d_Te,
     .                    d_Te_new,Sexp_Te,Simp1_Te_new,Simp2_Te_new,
     .                    v_Ti,a_Ti,a_Ti_new,b_Ti,c_Ti,c_Ti_new,d_Ti,
     .                    d_Ti_new,Sexp_Ti,Simp1_Ti_new,Simp2_Ti_new,
     .                    x_n,x_v,f_n,f_v,f_Te,f_Ti,f_n0,f_v0,f_T0,N,
     .                    mi,m0,B_K,Sext_n,Sext_v,Sext_Te,Sext_Ti,
     .                    Seir_n,Seir_v,Seir_Te,Seir_Ti,ION,REC,CX,
     .                    EXC,etai,kappae,kappai,exchange,Bsol)
      IMPLICIT NONE	
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::v_Te,v_Ti
      !velocity field
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::a_Te,a_Ti,a_Te_new,
     .                                             a_Ti_new,b_Te,b_Ti
      !coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::c_Te,c_Ti,c_Te_new,
     .                                             c_Ti_new,d_Te,d_Ti,
     .                                             d_Te_new,d_Ti_new
      !coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::Sexp_Te,Sexp_Ti,
     .                                             Simp1_Te_new,
     .                                             Simp1_Ti_new,
     .                                             Simp2_Te_new,
     .                                             Simp2_Ti_new
      !sources
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_n,x_v
      !grid
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::f_n,f_v,f_Te,f_Ti,f_n0,
     .                                          f_v0,f_T0
      !solved quantities
      INTEGER,INTENT(IN)::N
      !number of temperature points
      DOUBLE PRECISION,INTENT(IN)::mi,m0,B_K
      !constants
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Sext_n,Sext_v,Sext_Te,
     .                                          Sext_Ti,Seir_n,Seir_v,
     .                                          Seir_Te,Seir_Ti
      !external source of particles and energy
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::ION,REC,CX,EXC
      !collision rate coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::etai,kappae,kappai,
     .                                          exchange
      !transport coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Bsol
      !magnetic field
      DOUBLE PRECISION::skappae(N-1),skappai(N-1),s_v(N),s_v0(N),
     .                  sBsol(N-1)
      !staggered functions	
      DOUBLE PRECISION::dvdx_n(N),dBdx_n(N)
      !derivatives dv/dx, dB/dx
      INTEGER::i

      !staggered functions
      skappae=STAGGER1(kappae,N,x_n,x_v)
      skappai=STAGGER1(kappai,N,x_n,x_v)
      s_v=STAGGER2(f_v,N-1,x_n,x_v)
      s_v0=STAGGER2(f_v0,N-1,x_n,x_v)
      sBsol=STAGGER1(Bsol,N,x_n,x_v)

      !derivatives
      DO i=2,N-1
        dvdx_n(i)=(f_v(i)-f_v(i-1))/(x_v(i)-x_v(i-1))
        dBdx_n(i)=(sBsol(i)-sBsol(i-1))/(x_v(i)-x_v(i-1))
      END DO
      dvdx_n(1)=(f_v(2)-f_v(1))/(x_v(2)-x_v(1))
      dvdx_n(N)=(f_v(N-1)-f_v(N-2))/(x_v(N-1)-x_v(N-2))
      dBdx_n(1)=(sBsol(2)-sBsol(1))/(x_v(2)-x_v(1))
      dBdx_n(N)=(sBsol(N-1)-sBsol(N-2))/(x_v(N-1)-x_v(N-2))

      !ELECTRON EQUATION
	
      !velocity field
      v_Te=f_v

      !coefficients
      a_Te=1.0D0
      a_Te_new=1.0D0
      b_Te=1.0D0
      c_Te=-2.0D0*Bsol/(3.0D0*f_n)
      c_Te_new=c_Te
      d_Te=skappae/sBsol
      d_Te_new=d_Te
	
      !explicit source term
      Sexp_Te=2.0D0*(Sext_Te+Seir_Te)/(3.0D0*B_K*f_n) 
      Sexp_Te=Sexp_Te-2.0D0*f_n0*(ION*13.6D0+EXC)/3.0D0

      !implicit source terms
      Simp1_Te_new=-dvdx_n/3.0D0
      Simp1_Te_new=Simp1_Te_new-2.0D0*s_v*dBdx_n/(3.0D0*Bsol)
      Simp1_Te_new=Simp1_Te_new+(Sext_n+Seir_n)/f_n+exchange+f_n0*ION-
     .             f_n*REC
      !Simp11	
      Simp2_Te_new=-exchange !Simp21

      !ION EQUATION

      !velocity field
      v_Ti=f_v

      !coefficients
      a_Ti=1.0D0
      a_Ti_new=1.0D0
      b_Ti=1.0D0
      c_Ti=-2.0D0*Bsol/(3.0D0*f_n)
      c_Ti_new=c_Ti
      d_Ti=skappai/sBsol
      d_Ti_new=d_Ti

      !explicit source term
      Sexp_Ti=2.0D0*(Sext_Ti+Seir_Ti)/(3.0D0*B_K*f_n)
      Sexp_Ti=Sexp_Ti+f_n0*f_T0*(ION+CX)
      Sexp_Ti=Sexp_Ti+f_n0*m0*s_v0*s_v0*(ION+CX)/(3.0D0*B_K)
      Sexp_Ti=Sexp_Ti+m0*s_v*s_v*(-f_n*REC-f_n0*CX)/(3.0D0*B_K)
      Sexp_Ti=Sexp_Ti-2.0D0*s_v*mi*
     .        (f_n0*s_v0*ION+(s_v0-s_v)*f_n0*CX)/(3.0D0*B_K)
      Sexp_Ti=Sexp_Ti+s_v*s_v*mi*(f_n0*ION+f_n*REC)/(3.0D0*B_K)
      Sexp_Ti=Sexp_Ti+s_v*s_v*mi*(Sext_n+Seir_n)/(3.0D0*B_K*f_n)
      Sexp_Ti=Sexp_Ti-2.0D0*s_v*(Seir_v+Sext_v)/(3.0D0*B_K*f_n)
      Sexp_Ti=Sexp_Ti+2.0D0*etai/(3.0D0*B_K*f_n)*dvdx_n*dvdx_n
      Sexp_Ti=Sexp_Ti+2.0D0*etai*s_v*dBdx_n*dvdx_n/(3.0D0*B_K*f_n*Bsol)
      Sexp_Ti=Sexp_Ti+etai*s_v*s_v*dBdx_n*dBdx_n/
     .        (6.0D0*B_K*f_n*Bsol*Bsol)
	
      !implicit source terms
      Simp1_Ti_new=-dvdx_n/3.0D0
      Simp1_Ti_new=Simp1_Ti_new-2.0D0*s_v*dBdx_n/(3.0D0*Bsol)
      Simp1_Ti_new=Simp1_Ti_new+
     .             (Sext_n+Seir_n)/f_n+exchange+f_n0*(ION+CX)
      !Simp12
      Simp2_Ti_new=-exchange
      !Simp22

      END SUBROUTINE

      !--------------COEFFICIENTS FOR NEUTRAL CONTINUITY EQUATION--------------
      SUBROUTINE COEFFS_N0(v_n0,a_n0,a_n0_new,b_n0,Sexp_n0,Simp_n0_new,
     .                     f_n,f_Te,f_n0,f_v0,N,ION,REC)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::v_n0
      !velocity field	
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::a_n0,a_n0_new,b_n0
      !coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::Sexp_n0,Simp_n0_new
      !sources
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::f_n,f_Te,f_n0,f_v0
      !density, electron temperature and neutral density and velocity
      INTEGER,INTENT(IN)::N
      !number of density points	
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::ION,REC
      !rate coefficients
		
      !velocity field
      v_n0=f_v0

      !coefficients
      a_n0=1.0D0
      a_n0_new=a_n0
      b_n0=1.0D0
	
      !sources	
      Sexp_n0=f_n*f_n*REC 
      Simp_n0_new=f_n*ION 
		
      END SUBROUTINE

      !--------------COEFFICIENTS FOR NEUTRAL MOMENTUM EQUATION----------------
      SUBROUTINE COEFFS_V0(v_v0,a_v0,a_v0_new,b_v0,c_v0,c_v0_new,d_v0,
     .                     d_v0_new,Sexp_v0,Simp_v0_new,x_n,x_v,f_n,
     .                     f_v,f_Te,f_n0,f_v0,f_T0,N,mi,m0,B_K,REC,CX) 
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::v_v0
      !velocity field	
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::a_v0,a_v0_new,b_v0,
     .                                             c_v0,c_v0_new,d_v0,
     .                                             d_v0_new
      !coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::Sexp_v0,Simp_v0_new
      !sources
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_n,x_v
      !grid
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::f_n,f_v,f_Te,f_n0,f_v0,
     .                                          f_T0
      !solved quantities
      INTEGER,INTENT(IN)::N
      !number of velocity points
      DOUBLE PRECISION,INTENT(IN)::mi,m0,B_K
      !constants		
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::REC,CX
      !rate coefficients
      DOUBLE PRECISION::sREC(N),sCX(N),s_n(N),s_n0(N),s_v0(N+1)
      !staggered functions	
      DOUBLE PRECISION::derivation				
      INTEGER::i

      !staggered functions	
      sREC=STAGGER1(REC,N+1,x_n,x_v)
      sCX=STAGGER1(CX,N+1,x_n,x_v)
      s_n=STAGGER1(f_n,N+1,x_n,x_v)
      s_n0=STAGGER1(f_n0,N+1,x_n,x_v)
      s_v0=STAGGER2(f_v0,N,x_n,x_v)

      !velocity field (in midpoints of velocity grid)
      v_v0=0.5D0*s_v0

      !coefficients
      a_v0=1.0D0
      a_v0_new=a_v0
      b_v0=1.0D0	
      c_v0=0.0D0
      c_v0_new=c_v0	
      d_v0=0.0D0
      d_v0_new=d_v0

      !sources
      DO i=1,N 
        derivation=B_K*(f_n0(i+1)*f_T0(i+1)-
     .             f_n0(i)*f_T0(i))/(x_n(i+1)-x_n(i))
        Sexp_v0(i)=derivation              		              
      END DO
      Sexp_v0=-Sexp_v0/(m0*s_n0)

      Sexp_v0=Sexp_v0+s_n*s_n*f_v*sREC/s_n0+s_n*f_v*sCX
      Simp_v0_new=s_n*s_n*sREC/s_n0+s_n*sCX	
	
      END SUBROUTINE

      !----------------------------GRID GENERATION-----------------------------
      SUBROUTINE MAKEGRID(x_n,x_v,grid_parameter,L,n,NEUT)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::x_n,x_v
      !grid for density and velocity
      DOUBLE PRECISION,INTENT(IN)::grid_parameter
      !parameter of grid refinement
      DOUBLE PRECISION,INTENT(IN)::L
      !size of computational domain
      INTEGER,INTENT(IN)::n
      !number of velocity points
      INTEGER,INTENT(IN)::NEUT
      !switch for neutrals
      DOUBLE PRECISION::xdelta,gp
      INTEGER::i
	
      IF (NEUT.EQ.0) THEN !equidistant grid if no neutrals
        gp=1.0D0
      ELSE
        gp=grid_parameter
      ENDIF

      !construction of main grid (velocity grid)
      IF (n/2*2 .NE. n) THEN !for odd n
        xdelta=0.1D0
        x_v((n+1)/2)=0.0D0		
        DO i=1,(n+1)/2-1 
          x_v((n+1)/2+i)=x_v((n+1)/2+i-1)+xdelta
          x_v((n+1)/2-i)=x_v((n+1)/2-i+1)-xdelta
          xdelta=xdelta*gp
        END DO
        x_v=x_v*(0.5D0*L/x_v(n))
        x_v=x_v+x_v(n)
      ELSE                   !for even n
        xdelta=0.1D0
        x_v(n/2)=-0.5D0*xdelta !
        x_v(n/2+1)=0.5D0*xdelta 
        DO i=1,n/2-1
          xdelta=xdelta*gp
          x_v(n/2-i)=x_v(n/2-i+1)-xdelta
          x_v(n/2+1+i)=x_v(n/2+1+i-1)+xdelta			
        END DO
        x_v=x_v*(0.5D0*L/x_v(n))
        x_v=x_v+x_v(n)		
      ENDIF

      !construction of staggered grid (density grid in midpoints of the main grid)
      x_n(1)=0.0D0
      x_n(n+1)=L
      DO i=2,n
        x_n(i)=0.5D0*(x_v(i)+x_v(i-1))
      END DO

      x_n=x_n-0.5D0*L
      x_v=x_v-0.5D0*L

      DO i=1,n/2
        x_v(i)=-x_v(n+1-i)
        x_n(i)=-x_n(n+2-i)
      END DO

      END SUBROUTINE

      !---------------------CONVERSION TO STAGGERED MESH-----------------------
      FUNCTION STAGGER1(F,N,x_n,x_v)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::F        	!function for conversion
      INTEGER,INTENT(IN)::N					!number of points of function
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_n,x_v  	!grid		
      DOUBLE PRECISION,DIMENSION(N-1)::STAGGER1          	!function on staggered mesh
      INTEGER::i

      !from x_n to x_v (N to N-1)
      STAGGER1(1)=F(1)
      STAGGER1(N-1)=F(N)
      DO i=2,N-2       
        STAGGER1(i)=LIN_INTERP(x_v(i),x_n(i),x_n(i+1),F(i),F(i+1))
        !STAGGER1(i)=0.5D0*(F(i)+F(i+1))
      END DO

      END FUNCTION

      FUNCTION STAGGER2(F,N,x_n,x_v)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::F     	!function for conversion
      INTEGER,INTENT(IN)::N					!number of points of function	
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_n,x_v  	!grid	
      DOUBLE PRECISION,DIMENSION(N+1)::STAGGER2		!function on staggered mesh
      INTEGER::i

      !from x_v to x_n (N to N+1)
      STAGGER2(1)=F(1)
      STAGGER2(N+1)=F(N)
      DO i=2,N
        STAGGER2(i)=LIN_INTERP(x_n(i),x_v(i-1),x_v(i),F(i-1),F(i))
        !STAGGER2(i)=0.5D0*(F(i-1)+F(i))
      END DO

      END FUNCTION

      !--------------------------LINEAR INTERPOLATION--------------------------
      FUNCTION LIN_INTERP(x,x1,x2,f1,f2)
      IMPLICIT NONE
      DOUBLE PRECISION::LIN_INTERP !function at x
      DOUBLE PRECISION,INTENT(IN)::x,x1,x2,f1,f2		   

      LIN_INTERP=(f1*(x2-x)-f2*(x1-x))/(x2-x1)

      END FUNCTION

      !-------------------------------TIME STEP--------------------------------
      FUNCTION TIME_STEP(x_n,f_v)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::x_n,f_v  	!grid and velocity
      DOUBLE PRECISION::vmax,delta_x_min			!maximum velocity, minimum grid step
      DOUBLE PRECISION::TIME_STEP				!time step

      vmax=MAXVAL(ABS(f_v))
      delta_x_min=x_n(2)-x_n(1)
      TIME_STEP=delta_x_min/(5.0D0*vmax)
				
      END FUNCTION

      !------------------------------COLLISIONS--------------------------------
      SUBROUTINE ATOM_DATA(ION,REC,CX,EXC,f_n,f_Te,N,NEUT)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)::ION,REC,CX,EXC   	!rate coefficients
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::f_n,f_Te           	!density and electron temperature
      INTEGER,INTENT(IN)::N,NEUT					!number of points and switch for neutrals
      INTEGER::i

      IF (NEUT.EQ.1) THEN
        DO i=1,N
          ION(i)=IONIZATION(f_Te(i))
          !REC(i)=RECOMBINATION(f_Te(i))
          REC(i)=RECOMBINATION2(f_n(i),f_Te(i))
          CX(i)=CHARGE_EXCHANGE(f_Te(i))
          EXC(i)=EXCITATION(f_Te(i))
        END DO
      ENDIF

      IF (NEUT.EQ.0) THEN
        ION=0.0D0
        REC=0.0D0
        CX=0.0D0
        EXC=0.0D0
      ENDIF	

      END SUBROUTINE

      !------------------------------ATOMIC DATA-------------------------------
      FUNCTION IONIZATION(T)
      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN)::T		!temperature
      DOUBLE PRECISION::IONIZATION		!collision rate coefficient <sigma*v> [m3/s]
      DOUBLE PRECISION::TT,X,S
        
      TT=T
      IF (TT.LT.1.0D0) TT=1.0D0
      X=LOG10(TT)

      IF (TT.GE.20.0D0) THEN 
        S=-0.5151D0*X-2.563D0/X-5.231D0	
      ELSE
        S=-3.054D0*X-15.72D0*EXP(-X)+1.603D0*EXP(-X*X)
      ENDIF

      IONIZATION=10.0D0**(S-6.0D0)  

      END FUNCTION

      FUNCTION CHARGE_EXCHANGE(T)
      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN)::T   	!temperature
      DOUBLE PRECISION::CHARGE_EXCHANGE     	!<sigma*v> [m3/s]
      DOUBLE PRECISION::TT,S
 
      TT=T
      IF(TT.LT.1.0D0) TT=1.0D0
      S=-14.0D0+LOG10(TT)/3.0D0
      CHARGE_EXCHANGE=10.0D0**S		

      END FUNCTION

      FUNCTION RECOMBINATION(T)
      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN)::T        	!temperature
      DOUBLE PRECISION::RECOMBINATION		!<sigma*v> [m3/s]
      DOUBLE PRECISION::TT,X,PI,WYKL
      DOUBLE PRECISION,DIMENSION(10)::A
      INTEGER::I,LMAX

      A=(/-1.183D1,-5.534D0,8.714D0,-1.293D1,8.008D0,-1.772D0,0.0D0,
     .    0.0D0,0.0D0,0.0D0/)

      LMAX=5
      TT=T
      IF(TT.LT.1.0D0) TT=1.0D0
      X=LOG10(TT)/3.0D0
      PI=A(LMAX+1)
      DO I=LMAX+1,2,-1
        PI=A(I-1)+PI*X
      END DO
      WYKL=PI-6.0D0
      RECOMBINATION=10.0D0**WYKL

      END FUNCTION


      FUNCTION RECOMBINATION2(N,T) 
      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN)::N,T     	!density and temperature 
      DOUBLE PRECISION::RECOMBINATION2	!<sigma*v> [m3/s]
      DOUBLE PRECISION,DIMENSION(9,9)::A
      DOUBLE PRECISION::E,DNE,TT,REC,SUM,RN,RT,RREC,RDNE,RTE,RNJ,RTI
      INTEGER::I,J,I1,J1

      RECOMBINATION2=0.0
	
      END FUNCTION

      FUNCTION EXCITATION(T) !old version of excitation
      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN)::T       	!temperature
      DOUBLE PRECISION::EXCITATION          	!<sigma*v> [m3/s]
      DOUBLE PRECISION::Y,TT
    
      TT=T
      IF(TT.LT.1.0D0) TT=1.0D0 
      Y=10.2D0/TT

      EXCITATION=49.0D-14/(0.28D0+Y)*EXP(-Y)*SQRT(Y*(1.0D0+Y))   			
	
      END FUNCTION

      !------------------------------OUTPUT FINAL------------------------------
      SUBROUTINE OUTPUT_FINAL(f_n,f_v,f_Te,f_Ti,f_n0,f_v0,x_n,x_v,
     .                        Sext_n,Sext_v,Sext_Te,Sext_Ti,Seir_n,
     .                        Seir_v,Seir_Te,Seir_Ti,Bsol,n)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::f_n,f_v,f_Te,f_Ti,f_n0,
     .                                          f_v0,x_n,x_v
      !quantities and grid      
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Sext_n,Sext_v,Sext_Te,
     .                                          Sext_Ti,Seir_n,Seir_v,
     .                                          Seir_Te,Seir_Ti
      !sources
      DOUBLE PRECISION,DIMENSION(:),INTENT(IN)::Bsol
      !magnetic field
      INTEGER,INTENT(IN)::n
      !number of grid points
      DOUBLE PRECISION::s_v(n+1),s_v0(n+1)	
      INTEGER::i

      !staggered quantities
      s_v=STAGGER2(f_v,n,x_n,x_v)
      s_v0=STAGGER2(f_v0,n,x_n,x_v)

      !output.dat
      OPEN(UNIT=11,FILE='output.dat',STATUS='REPLACE',ACTION='WRITE')
      WRITE(11,fmt="(1I)") n
      DO i=1,n+1
        WRITE(11,fmt="(16E)") x_n(i),f_n(i),s_v(i),f_Te(i),f_Ti(i),
     .                        f_n0(i),s_v0(i),Bsol(i),Sext_n(i),
     .                        Sext_v(i),Sext_Te(i),Sext_Ti(i),Seir_n(i),
     .                        Seir_v(i),Seir_Te(i),Seir_Ti(i)
      END DO
      CLOSE(UNIT=11)	

      END SUBROUTINE

      END MODULE SOLF1D

      !------------------------------END OF CODE-------------------------------




