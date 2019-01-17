      SUBROUTINE Solve_Potential(ion)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none

c    This routine is used to solve Ohm's law and the current conservation equation in OSM.
c    The two first-order PDEs are combined into one 2nd order PDE, and solved as a standard
c    1-D convection-diffusion equation using a Hybrid discretisation method.
c
c    Revision history:
c     - 1.0 (5th January 2010) First go at writing the solver. There was a sign problem
c       in the equation, fixed by multiplying De on the inlet boundary by -1.
c
c     - 1.1 (26th July 2010)   Re-derivation of the numerical scheme resolved the above problem.
c       it was also found that the radial electric field using solutions from the earlier
c       version of the solver were unrealistically large. This has now been resolved.
c       The solver has also been cross-checked against a TVD scheme outlined in Versteeg,
c       and shown to be in good agreement (better than 1% across the solution).
c
c     - 1.2 (2nd August 2010)  Had a go at re-deriving the Hybrid solver again, as testing
c       was showing that it wasn't switching between FOU and 2nd order central differencing
c       properly. This has now been fixed, and the electrostatic potential solutions
c       from the code now agree very well with a modified MATLAB code (check_efield3)
c       which has been successfully benchmarked.
c
c     - 1.3 (17th August 2010) Found another bug, the code now agrees perfectly, not very
c       well, but perfectly, with MATLAB.

      INTEGER er, ic, i, ion, ict

      REAL*8  nL,nR,gL,tL,tR,pL,pR,m_e,dp,dm,
     .        beta,tmp,etaL,etaR,phi_L,phi_R

      REAL*8  eta(icmax),etap(icmax),nep(icmax),tep(icmax),
     .        tepp(icmax),teppc(icmax),p(icmax),pp(icmax),
     .        ppp(icmax),gradbovb(icmax),delta_s(icmax),
     .        Sj(icmax),LHS_TDMA(3,icmax),RHS(icmax),
     .        Fe(icmax),Fw(icmax),De(icmax),Dw(icmax),
     .        big_F(icmax),big_D(icmax),LHS_tmp(3)

c     Declare some constants
c     TODO: These should really take values from OSM...
      m_e   = 9.10938188D-31

c     Start to fill up these arrays...

      DO i = 1,icmax
        Sj(i)       = 0.0D0 !cross-field current source
        eta(i)      = 1.0D0/(3.6D7*(te(i)/1000.0D0)**1.5D0)
        p(i)        = ne(i)*te(i)
        delta_s(i)  = cell(i)%sbnd(2)-cell(i)%sbnd(1)
        epot(i)      = 0.0D0
      ENDDO

c     Calculate the first derivatives

      DO i = 2,(icmax-1)
        dp   = SFOR(i+1)-SFOR(i)
        dm   = SFOR(i)-SFOR(i-1)
        beta = dp/dm

        nep(i) = (ne(i+1)-(1.0D0-beta**2.0D0)*ne(i)
     .          -(beta**2.0D0)*ne(i-1))/(dp*(1.0D0+beta))
        pp(i)  = (p(i+1)-(1.0D0-beta**2.0D0)*p(i)
     .          -(beta**2.0D0)*p(i-1))/(dp*(1.0D0+beta))
        tep(i) = (te(i+1)-(1.0D0-beta**2.0D0)*te(i)
     .          -(beta**2.0D0)*te(i-1))/(dp*(1.0D0+beta))
        etap(i)= (eta(i+1)-(1.0D0-beta**2.0D0)*eta(i)
     .          -(beta**2.0D0)*eta(i-1))/(dp*(1.0D0+beta))
        gradbovb(i) = (bfield(i+1)-(1.0D0-beta**2.0D0)*bfield(i)
     .          -(beta**2.0D0)*bfield(i-1))/(dp*(1.0D0+beta))
        gradbovb(i) = gradbovb(i)/bfield(i)
      ENDDO

c     Calculate the second derivatives

      DO i = 2,(icmax-1)
        dp   = SFOR(i+1)-SFOR(i)
        dm   = SFOR(i)-SFOR(i-1)
        beta = dp/dm

        tepp(i) = (2.0D0*(te(i+1)-(1.0D0+beta)*te(i)
     .             +beta*te(i-1)))/(dp*dm*(1.0D0+beta))
        ppp(i)  = (2.0D0*(p(i+1)-(1.0D0+beta)*p(i)
     .             +beta*p(i-1)))/(dp*dm*(1.0D0+beta))
      ENDDO

      nL     = ne(0)
      nR     = ne(icmax+1)
      tL     = te(0)
      tR     = te(icmax+1)
      pL     = nL*tL
      pR     = nR*tR
      etaL   = 1.0D0/(3.6D7*(tL/1000.0D0)**1.5D0)
      etaR   = 1.0D0/(3.6D7*(tR/1000.0D0)**1.5D0)

c     Calculate the first derivative at the edges

c    Start of the domain

      dp   = SFOR(2)-SFOR(1)
      dm   = SFOR(1)
      beta = dp/dm

      nep(1) = (ne(2)-(1.0D0-beta**2.0D0)*ne(1)-(beta**2.0D0)*nL)/
     .         (dp*(1.0D0+beta))

      pp(1) = (p(2)-(1.0D0-beta**2.0D0)*p(1)-(beta**2.0D0)*pL)/
     .         (dp*(1.0D0+beta))

      tep(1) = (te(2)-(1.0D0-beta**2.0D0)*te(1)-(beta**2.0D0)*tL)/
     .         (dp*(1.0D0+beta))

      etap(1) = (eta(2)-(1.0D0-beta**2.0D0)*eta(1)-(beta**2.0D0)*etaL)/
     .         (dp*(1.0D0+beta))
     
c      gradbovb(1) = (bfield(2)-(1.0D0-beta**2.0D0)*bfield(1)-
c     .              (beta**2.0D0)*bfield(1))/
c     .              (dp*(1.0D0+beta))
     
c      gradbovb(1) = gradbovb(1)/bfield(1)

c     Use second cell value in grid
      gradbovb(1) = gradbovb(2)

      tepp(1) = 2.0D0*(te(2)-(1.0D0+beta)*te(1)+beta*tL)/
     .                (dp*dm*(1.0D0+beta))

      ppp(1)  = 2.0D0*(p(2)-(1.0D0+beta)*p(1)+beta*pL)/
     .                (dp*dm*(1.0D0+beta))

c     End of the domain

      dp = cell(icmax)%sbnd(2)-SFOR(icmax)
      dm = SFOR(icmax)-SFOR(icmax-1)
      beta = dp/dm

      nep(icmax) = (nR-(1.0D0-beta**2.0D0)*ne(icmax)-
     .             (beta**2.0D0)*ne(icmax-1))/
     .             (dp*(1.0D0+beta))

      pp(icmax)  = (pR-(1.0D0-beta**2.0D0)*p(icmax)-
     .             (beta**2.0D0)*p(icmax-1))/
     .             (dp*(1.0D0+beta))

      tep(icmax) = (tR-(1.0D0-beta**2.0D0)*te(icmax)-
     .             (beta**2.0D0)*te(icmax-1))/
     .             (dp*(1.0D0+beta))

      etap(icmax) = (etaR-(1.0D0-beta**2.0D0)*eta(icmax)-
     .             (beta**2.0D0)*eta(icmax-1))/
     .             (dp*(1.0D0+beta))
     
c      gradbovb(icmax) = (bfield(icmax)-
c     .             (1.0D0-beta**2.0D0)*bfield(icmax)-
c     .             (beta**2.0D0)*bfield(icmax-1))/
c     .             (dp*(1.0D0+beta))
c      gradbovb(icmax) = gradbovb(icmax)/bfield(icmax)
      
c     Use previous cell value...
      gradbovb(icmax) = gradbovb(icmax-1)

      tepp(icmax) = 2.0D0*(tR-(1.0D0+beta)*te(icmax)+
     .                     beta*te(icmax-1))/
     .                    (dp*dm*(1.0D0+beta))

      ppp(icmax)  = 2.0D0*(pR-(1.0D0+beta)*p(icmax)+
     .                     beta*p(icmax-1))/
     .                    (dp*dm*(1.0D0+beta))

      DO i = 1,icmax
        big_F(i) = -1.0D0*gradbovb(i)-(1.0D0/eta(i))*etap(i)
        big_D(i) = 1.0D0
      ENDDO

c     Set boundary conditions - floating walls, see
c     Stageby p.79
c     TODO: THIS SHOULD BE AN OPTION IN OSM INPUT FILE!

      phi_L = 0.5D0*tL*LOG(4.0D0*3.141592D0*(m_e/mi(ion)))
      phi_R = 0.5D0*tR*LOG(4.0D0*3.141592D0*(m_e/mi(ion)))

c     Calculate contents of Fe,Fw,De,Dw

      Fw(1) = big_F(1)-(delta_s(1)/(delta_s(2)+delta_s(1)))
     .         *(big_F(2)-big_F(1))

      Fe(1) = big_F(1)+(delta_s(1)/(delta_s(2)+delta_s(1)))
     .         *(big_F(2)-big_F(1))

      Dw(1) = big_D(1)-(delta_s(1)/(delta_s(2)+delta_s(1)))
     .         *(big_D(2)-big_D(1))

      Dw(1) = (2*Dw(1))/delta_s(1)

      De(1) = big_D(1)+(delta_s(1)/(delta_s(2)+delta_s(1)))
     .         *(big_D(2)-big_D(1))

      De(1) = (2*De(1))/(delta_s(1)+delta_s(2))

      DO i = 2,(icmax-1)
        dp = delta_s(i-1)/(delta_s(i-1)+delta_s(i))
        dm = delta_s(i)/(delta_s(i)+delta_s(i+1))

        Fw(i) = (1.0D0-dp)*big_F(i-1)+dp*big_F(i)
        Fe(i) = (1.0D0-dm)*big_F(i)+dm*big_F(i+1)

        Dw(i) = (1.0D0-dp)*big_D(i-1)+dp*big_D(i)
        Dw(i) = (2.0D0*Dw(i))/(delta_s(i-1)+delta_s(i))
        De(i) = (1.0D0-dm)*big_D(i)+dm*big_D(i+1)
        De(i) = (2.0D0*De(i))/(delta_s(i)+delta_s(i+1))
      ENDDO

      Fw(icmax) = big_F(icmax-1)+(delta_s(icmax-1)
     .           /(delta_s(icmax)+delta_s(icmax-1)))
     .           *(big_F(icmax)-big_F(icmax-1))

      Fe(icmax) = big_F(icmax-1)+((delta_s(icmax-1)
     .            +2.0D0*delta_s(icmax))/
     .            (delta_s(icmax)+delta_s(icmax-1)))
     .            *(big_F(icmax)-big_F(icmax-1))

      Dw(icmax) = big_D(icmax-1)+(delta_s(icmax-1)
     .            /(delta_s(icmax)+delta_s(icmax-1)))
     .            *(big_D(icmax)-big_D(icmax-1))

      Dw(icmax) = (2.0D0*Dw(icmax))/
     .            (delta_s(icmax-1)+delta_s(icmax))

      De(icmax) = big_D(icmax-1)+
     .            ((delta_s(icmax-1)+2.0D0*delta_s(icmax))/
     .            (delta_s(icmax)+delta_s(icmax-1)))
     .            *(big_D(icmax)-big_D(icmax-1))

      De(icmax) = (2.0D0*De(icmax))/(delta_s(icmax))

c     Calculate the value of RHS (take a deep breath...)

      DO i = 1,icmax
        RHS(i) = (1.0D0/ne(i))*ppp(i)
     .          -(1.0D0/(ne(i)**2.0D0))*pp(i)*nep(i)
     .          -(1.0D0/ne(i))*gradbovb(i)*pp(i)
     .          -(1.0D0/(ne(i)*eta(i)))*pp(i)*etap(i)
     .          +0.71D0*tepp(i)
     .          -0.71D0*tep(i)*gradbovb(i)
     .          -(0.71D0/eta(i))*tep(i)*etap(i)
     .          -eta(i)*Sj(i)

        RHS(i) = RHS(i)*(1.0D0)*delta_s(i)
      ENDDO

c     Calculate the LHS (using the Hybrid scheme)

c     Inlet boundary
      IF (big_F(1).LT.0.0D0) THEN
        LHS_TDMA(1,1) = 0.0D0
        LHS_TDMA(2,1) = -1.0D0*Fw(1)-De(1)-Dw(1)
        LHS_TDMA(1,1) = De(1)+Fe(1)
        RHS(1)        = RHS(1)-phi_L*Dw(1)
      ELSE
        LHS_TDMA(1,1) = 0.0D0
        LHS_TDMA(2,1) = (Fe(1)-De(1)-Dw(1))
        LHS_TDMA(3,1) = De(1)
        RHS(1)        = RHS(1)-phi_L*(Dw(1)-Fw(1))
      ENDIF

      DO i = 2,(icmax-1)
        dm = delta_s(i-1)/(delta_s(i-1)+delta_s(i))
        dp = delta_s(i)/(delta_s(i)+delta_s(i+1))

c       Left cell contribution (i)
        LHS_tmp(1)    = -1.0D0*Fw(i)
        LHS_tmp(2)    = Dw(i)-Fw(i)*(1.0D0-dm)
        LHS_tmp(3)    = 0.0D0
        IF(DABS(Fw(i)/Dw(i)).LT.2.0D0)THEN
          LHS_TDMA(1,i) = LHS_tmp(2)
        ELSE
          LHS_TDMA(1,i) = Dw(i)+MINVAL([LHS_tmp(1),LHS_tmp(3)])
        ENDIF

c       Right cell contribution (i)
        LHS_tmp(1)    = 1.0D0*Fe(i)
        LHS_tmp(2)    = De(i)+dp*Fe(i)!dp*Fe(i)
        LHS_tmp(3)    = 0.0D0
        IF(DABS(Fe(i)/De(i)).LT.2.0D0)THEN
          LHS_TDMA(3,i) = LHS_tmp(2)
        ELSE
          LHS_TDMA(3,i) = De(i)+MINVAL([LHS_tmp(1),LHS_tmp(3)])
        ENDIF

c       Central cell contribution (i)
        LHS_TDMA(2,i) = (Fe(i)-Fw(i))
     .                  -LHS_TDMA(3,i)
     .                  -LHS_TDMA(1,i)

      ENDDO

c     Outlet boundary
      IF (big_F(icmax).LT.0.0D0) THEN
        LHS_TDMA(2,icmax) = -1.0D0*Fw(icmax)-De(icmax)-Dw(icmax)
        LHS_TDMA(1,icmax) = Dw(icmax)
        RHS(icmax) = RHS(icmax)-phi_R*(De(icmax)+Fe(icmax))
      ELSE
        LHS_TDMA(2,icmax) = Fe(icmax)-De(icmax)-Dw(icmax)
        LHS_TDMA(1,icmax) = Dw(icmax)-Fw(icmax)
        RHS(icmax) = RHS(icmax)-phi_R*(De(icmax)+Fe(icmax))
      ENDIF

c     Solve for phi using the TDMA
      CALL tridag(LHS_TDMA(1,:),LHS_TDMA(2,:),LHS_TDMA(3,:),
     .            RHS,epot(1:icmax),icmax)

      epot(0) = phi_L
      epot(icmax+1) = phi_R

c     Calculate the electrostatic field
      DO i = 2,(icmax-1)
        dp = SFOR(i+1)-SFOR(i)
        dm = SFOR(i)-SFOR(i-1)
        beta = dp/dm

        efield(i) = (epot(i+1)-(1-beta**2.0D0)*epot(i)
     .               -(beta**2.0D0)*epot(i-1))
     .               /(dp*(1.0D0+beta))
      ENDDO

      efield(1)     = (epot(1)-phi_L)/
     .                (SFOR(1)-cell(1)%sbnd(1))
      efield(icmax) = (phi_R-epot(icmax))/
     .                (cell(icmax)%sbnd(2)-SFOR(icmax))

      efield(:) = -1.0D0*efield(:)

      epot(0)   = phi_L
      epot(icmax+1) = phi_R

      END

      SUBROUTINE tridag(a,b,c,r,u,n)
      INTEGER n,NMAX
      REAL*8 a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=500)
c This routine solves a tridiagonal set of linear equations, the algorithm was
c lifted from Numerical Recipes in Fortran, p. 43.  It was modified slightly
c changing real to real*8 for compatibility with OSM.

c Solves for a vector u(1:n) of length n the tridiagonal linear set given by equation (2.4.1).
c a(1:n), b(1:n), c(1:n), and r(1:n) are input vectors and are not modified.
c Parameter: NMAX is the maximum expected value of n.
      INTEGER j
      REAL*8 bet,gam(NMAX)  ! One vector of workspace, gam is needed.
      if(b(1).eq.0.)pause 'tridag: rewrite equations'
c If this happens then you should rewrite your equations as a set of order N − 1, with u2
c trivially eliminated.
      bet=b(1)
      u(1)=r(1)/bet
      do j=2,n ! Decomposition and forward substitution.
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.)pause 'tridag failed' ! Algorithm fails; see below.
        u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
      do j=n-1,1,-1 ! Backsubstitution.
        u(j)=u(j)-gam(j+1)*u(j+1)
      enddo
      return
      END