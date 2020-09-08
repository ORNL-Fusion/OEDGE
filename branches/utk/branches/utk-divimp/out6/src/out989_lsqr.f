*     From arpa!sol-michael.stanford.edu!mike 5 May 89 23:53:00 PDT
*     Original source downloaded from www.netlib.org (search for LSQR) - SL 29.09.05
      SUBROUTINE LSQR  ( M, N, APROD, DAMP,
     $                   LENIW, LENRW, IW, RW,
     $                   U, V, W, X, SE,
     $                   ATOL, BTOL, CONLIM, ITNLIM, NOUT,
     $                   ISTOP, ITN, ANORM, ACOND, RNORM, ARNORM, XNORM)

      implicit none
      EXTERNAL           APROD
      INTEGER            M, N, LENIW, LENRW, ITNLIM, NOUT, ISTOP, ITN
      INTEGER            IW(LENIW)
      DOUBLE PRECISION   RW(LENRW), U(M), V(N), W(N), X(N), SE(N),
     $                   ATOL, BTOL, CONLIM, DAMP,
     $                   ANORM, ACOND, RNORM, ARNORM, XNORM,
c slmod begin
c ----------------------------------------------------------------------
     $                   OLDX(N)
c ----------------------------------------------------------------------
c slmod end
*-----------------------------------------------------------------------
*
*     LSQR  finds a solution x to the following problems:
*
*     1. Unsymmetric equations --    solve  A*x = b
*
*     2. Linear least squares  --    solve  A*x = b
*                                    in the least-squares sense
*
*     3. Damped least squares  --    solve  (   A    )*x = ( b )
*                                           ( damp*I )     ( 0 )
*                                    in the least-squares sense
*
*     where A is a matrix with m rows and n columns, b is an
*     m-vector, and damp is a scalar.  (All quantities are real.)
*     The matrix A is intended to be large and sparse.  It is accessed
*     by means of subroutine calls of the form
*
*                CALL APROD ( mode, m, n, x, y, LENIW, LENRW, IW, RW )
*
*     which must perform the following functions:
*
*                If MODE = 1, compute  y = y + A*x.
*                If MODE = 2, compute  x = x + A(transpose)*y.
*
*     The vectors x and y are input parameters in both cases.
*     If  mode = 1,  y should be altered without changing x.
*     If  mode = 2,  x should be altered without changing y.
*     The parameters LENIW, LENRW, IW, RW may be used for workspace
*     as described below.
*
*     The rhs vector b is input via U, and subsequently overwritten.
*
*
*     Note:  LSQR uses an iterative method to approximate the solution.
*     The number of iterations required to reach a certain accuracy
*     depends strongly on the scaling of the problem.  Poor scaling of
*     the rows or columns of A should therefore be avoided where
*     possible.
*
*     For example, in problem 1 the solution is unaltered by
*     row-scaling.  If a row of A is very small or large compared to
*     the other rows of A, the corresponding row of ( A  b ) should be
*     scaled up or down.
*
*     In problems 1 and 2, the solution x is easily recovered
*     following column-scaling.  Unless better information is known,
*     the nonzero columns of A should be scaled so that they all have
*     the same Euclidean norm (e.g., 1.0).
*
*     In problem 3, there is no freedom to re-scale if damp is
*     nonzero.  However, the value of damp should be assigned only
*     after attention has been paid to the scaling of A.
*
*     The parameter damp is intended to help regularize
*     ill-conditioned systems, by preventing the true solution from
*     being very large.  Another aid to regularization is provided by
*     the parameter ACOND, which may be used to terminate iterations
*     before the computed solution becomes very large.
*
*
*     Notation
*     --------
*
*     The following quantities are used in discussing the subroutine
*     parameters:
*
*     Abar   =  (   A    ),          bbar  =  ( b )
*               ( damp*I )                    ( 0 )
*
*     r      =  b  -  A*x,           rbar  =  bbar  -  Abar*x
*
*     rnorm  =  sqrt( norm(r)**2  +  damp**2 * norm(x)**2 )
*            =  norm( rbar )
*
*     RELPR  =  the relative precision of floating-point arithmetic
*               on the machine being used.  For example, on the IBM 370,
*               RELPR is about 1.0E-6 and 1.0D-16 in single and double
*               precision respectively.
*
*     LSQR  minimizes the function rnorm with respect to x.
*
*
*     Parameters
*     ----------
*
*     M       input      m, the number of rows in A.
*
*     N       input      n, the number of columns in A.
*
*     APROD   external   See above.
*
*     DAMP    input      The damping parameter for problem 3 above.
*                        (DAMP should be 0.0 for problems 1 and 2.)
*                        If the system A*x = b is incompatible, values
*                        of DAMP in the range 0 to sqrt(RELPR)*norm(A)
*                        will probably have a negligible effect.
*                        Larger values of DAMP will tend to decrease
*                        the norm of x and reduce the number of
*                        iterations required by LSQR.
*
*                        The work per iteration and the storage needed
*                        by LSQR are the same for all values of DAMP.
*
*     LENIW   input      The length of the workspace array IW.
*     LENRW   input      The length of the workspace array RW.
*     IW      workspace  An integer array of length LENIW.
*     RW      workspace  A real array of length LENRW.
*
*             Note:  LSQR  does not explicitly use the previous four
*             parameters, but passes them to subroutine APROD for
*             possible use as workspace.  If APROD does not need
*             IW or RW, the values LENIW = 1 or LENRW = 1 should
*             be used, and the actual parameters corresponding to
*             IW or RW  may be any convenient array of suitable type.
*
*     U(M)    input      The rhs vector b.  Beware that U is
*                        over-written by LSQR.
*
*     V(N)    workspace
*     W(N)    workspace
*
*     X(N)    output     Returns the computed solution x.
*
*     SE(N)   output     Returns standard error estimates for the
*                        components of X.  For each i, SE(i) is set
*                        to the value  rnorm * sqrt( sigma(i,i) / T ),
*                        where sigma(i,i) is an estimate of the i-th
*                        diagonal of the inverse of Abar(transpose)*Abar
*                        and  T = 1      if  m .le. n,
*                             T = m - n  if  m .gt. n  and  damp = 0,
*                             T = m      if  damp .ne. 0.
*
*     ATOL    input      An estimate of the relative error in the data
*                        defining the matrix A.  For example,
*                        if A is accurate to about 6 digits, set
*                        ATOL = 1.0E-6 .
*
*     BTOL    input      An extimate of the relative error in the data
*                        defining the rhs vector b.  For example,
*                        if b is accurate to about 6 digits, set
*                        BTOL = 1.0E-6 .
*
*     CONLIM  input      An upper limit on cond(Abar), the apparent
*                        condition number of the matrix Abar.
*                        Iterations will be terminated if a computed
*                        estimate of cond(Abar) exceeds CONLIM.
*                        This is intended to prevent certain small or
*                        zero singular values of A or Abar from
*                        coming into effect and causing unwanted growth
*                        in the computed solution.
*
*                        CONLIM and DAMP may be used separately or
*                        together to regularize ill-conditioned systems.
*
*                        Normally, CONLIM should be in the range
*                        1000 to 1/RELPR.
*                        Suggested value:
*                        CONLIM = 1/(100*RELPR)  for compatible systems,
*                        CONLIM = 1/(10*sqrt(RELPR)) for least squares.
*
*             Note:  If the user is not concerned about the parameters
*             ATOL, BTOL and CONLIM, any or all of them may be set
*             to zero.  The effect will be the same as the values
*             RELPR, RELPR and 1/RELPR respectively.
*
*     ITNLIM  input      An upper limit on the number of iterations.
*                        Suggested value:
*                        ITNLIM = n/2   for well-conditioned systems
*                                       with clustered singular values,
*                        ITNLIM = 4*n   otherwise.
*
*     NOUT    input      File number for printed output.  If positive,
*                        a summary will be printed on file NOUT.
*
*     ISTOP   output     An integer giving the reason for termination:
*
*                0       x = 0  is the exact solution.
*                        No iterations were performed.
*
*                1       The equations A*x = b are probably
*                        compatible.  Norm(A*x - b) is sufficiently
*                        small, given the values of ATOL and BTOL.
*
*                2       The system A*x = b is probably not
*                        compatible.  A least-squares solution has
*                        been obtained that is sufficiently accurate,
*                        given the value of ATOL.
*
*                3       An estimate of cond(Abar) has exceeded
*                        CONLIM.  The system A*x = b appears to be
*                        ill-conditioned.  Otherwise, there could be an
*                        error in subroutine APROD.
*
*                4       The equations A*x = b are probably
*                        compatible.  Norm(A*x - b) is as small as
*                        seems reasonable on this machine.
*
*                5       The system A*x = b is probably not
*                        compatible.  A least-squares solution has
*                        been obtained that is as accurate as seems
*                        reasonable on this machine.
*
*                6       Cond(Abar) seems to be so large that there is
*                        no point in doing further iterations,
*                        given the precision of this machine.
*                        There could be an error in subroutine APROD.
*
*                7       The iteration limit ITNLIM was reached.
*
*     ITN     output     The number of iterations performed.
*
*     ANORM   output     An estimate of the Frobenius norm of  Abar.
*                        This is the square-root of the sum of squares
*                        of the elements of Abar.
*                        If DAMP is small and if the columns of A
*                        have all been scaled to have length 1.0,
*                        ANORM should increase to roughly sqrt(n).
*                        A radically different value for ANORM may
*                        indicate an error in subroutine APROD (there
*                        may be an inconsistency between modes 1 and 2).
*
*     ACOND   output     An estimate of cond(Abar), the condition
*                        number of Abar.  A very high value of ACOND
*                        may again indicate an error in APROD.
*
*     RNORM   output     An estimate of the final value of norm(rbar),
*                        the function being minimized (see notation
*                        above).  This will be small if A*x = b has
*                        a solution.
*
*     ARNORM  output     An estimate of the final value of
*                        norm( Abar(transpose)*rbar ), the norm of
*                        the residual for the usual normal equations.
*                        This should be small in all cases.  (ARNORM
*                        will often be smaller than the true value
*                        computed from the output vector X.)
*
*     XNORM   output     An estimate of the norm of the final
*                        solution vector X.
*
*
*     Subroutines and functions used
*     ------------------------------
*
*     USER               APROD
*     BLAS               DCOPY, DNRM2, DSCAL (see Lawson et al. below)
*
*
*     Precision
*     ---------
*
*     The number of iterations required by LSQR will usually decrease
*     if the computation is performed in higher precision.  To convert
*     LSQR between single and double precision, change the words
*                        DOUBLE PRECISION
*                        DCOPY, DNRM2, DSCAL
*     to the appropriate FORTRAN and BLAS equivalents.
*     Also change 'D+' or 'E+' in the PARAMETER statement.
*
*
*     References
*     ----------
*
*     C.C. Paige and M.A. Saunders,  LSQR: An algorithm for sparse
*          linear equations and sparse least squares,
*          ACM Transactions on Mathematical Software 8, 1 (March 1982),
*          pp. 43-71.
*
*     C.C. Paige and M.A. Saunders,  Algorithm 583, LSQR: Sparse
*          linear equations and least-squares problems,
*          ACM Transactions on Mathematical Software 8, 2 (June 1982),
*          pp. 195-209.
*
*     C.L. Lawson, R.J. Hanson, D.R. Kincaid and F.T. Krogh,
*          Basic linear algebra subprograms for Fortran usage,
*          ACM Transactions on Mathematical Software 5, 3 (Sept 1979),
*          pp. 308-323 and 324-325.
*-----------------------------------------------------------------------
*
*
*     LSQR development:
*     22 Feb 1982: LSQR sent to ACM TOMS to become Algorithm 583.
*     15 Sep 1985: Final F66 version.  LSQR sent to "misc" in netlib.
*     13 Oct 1987: Bug (Robert Davies, DSIR).  Have to delete
*                     IF ( (ONE + DABS(T)) .LE. ONE ) GO TO 200
*                  from loop 200.  The test was an attempt to reduce
*                  underflows, but caused W(I) not to be updated.
*     17 Mar 1989: First F77 version.
*     04 May 1989: Bug (David Gay, AT&T).  When the second BETA is zero,
*                  RNORM = 0 and
*                  TEST2 = ARNORM / (ANORM * RNORM) overflows.
*                  Fixed by testing for RNORM = 0.
*     05 May 1989: Sent to "misc" in netlib.
*
*     Michael A. Saunders            (na.saunders @ NA-net.stanford.edu)
*     Department of Operations Research
*     Stanford University
*     Stanford, CA 94305-4022.
*-----------------------------------------------------------------------

*     Intrinsics and local variables

      INTRINSIC          ABS, MOD, SQRT
      INTEGER            I, NCONV, NSTOP
      DOUBLE PRECISION   DNRM2
      DOUBLE PRECISION   ALFA, BBNORM, BETA, BNORM,
     $                   CS, CS1, CS2, CTOL, DAMPSQ, DDNORM, DELTA,
     $                   GAMMA, GAMBAR, PHI, PHIBAR, PSI,
     $                   RES1, RES2, RHO, RHOBAR, RHBAR1, RHBAR2,
     $                   RHS, RTOL, SN, SN1, SN2,
     $                   T, TAU, TEST1, TEST2, TEST3,
     $                   THETA, T1, T2, T3, XXNORM, Z, ZBAR

      DOUBLE PRECISION   ZERO,           ONE
      PARAMETER        ( ZERO = 0.0D+0,  ONE = 1.0D+0 )

      CHARACTER*16       ENTER, EXIT
      CHARACTER*60       MSG(0:8)

      DATA               ENTER /' Enter LSQR.    '/,
     $                   EXIT  /' Exit  LSQR.    '/

      DATA               MSG
     $ / 'The exact solution is  X = 0',
     $   'Ax - b is small enough, given ATOL, BTOL',
     $   'The least-squares solution is good enough, given ATOL',
     $   'The estimate of cond(Abar) has exceeded CONLIM',
     $   'Ax - b is small enough for this machine',
     $   'The least-squares solution is good enough for this machine',
     $   'Cond(Abar) seems to be too large for this machine',
     $   'The iteration limit has been reached',
     $   'Negative x value detected' /
c slmod begin
c ----------------------------------------------------------------------
c      INTEGER initialize
c      DATA    initialize / 0 /

       LOGICAL STOP8OKAY

c      SAVE

c      IF (initialize.EQ.1) GOTO 90

c      WRITE(6,*) 'INITIALIZING LSQR'

c      initialize = 1

c      WRITE(0,*) 'ITNLIM=',itnlim
      STOP8OKAY=.FALSE.
      IF (ITNLIM.LT.0) THEN
        ITNLIM = -ITNLIM
        WRITE(0,*) 'WATCHING FOR NEGATIVE VALUES'
        STOP8OKAY = .TRUE.
      ENDIF
c      WRITE(0,*) 'STOP8OKAY',stop8okay

c ----------------------------------------------------------------------
c slmod end
*-----------------------------------------------------------------------


*     Initialize.

      IF (NOUT .GT. 0)
     $   WRITE(NOUT, 1000) ENTER, M, N, DAMP, ATOL, CONLIM, BTOL, ITNLIM
      ITN    =   0
      ISTOP  =   0
      NSTOP  =   0
      CTOL   =   ZERO
      IF (CONLIM .GT. ZERO) CTOL = ONE / CONLIM
      ANORM  =   ZERO
      ACOND  =   ZERO
      BBNORM =   ZERO
      DAMPSQ =   DAMP**2
      DDNORM =   ZERO
      RES2   =   ZERO
      XNORM  =   ZERO
      XXNORM =   ZERO
      CS2    = - ONE
      SN2    =   ZERO
      Z      =   ZERO

      DO 10  I = 1, N
         V(I)  =  ZERO
         X(I)  =  ZERO
        SE(I)  =  ZERO
   10 CONTINUE

*     Set up the first vectors U and V for the bidiagonalization.
*     These satisfy  BETA*U = b,  ALFA*V = A(transpose)*U.

      ALFA   =   ZERO
      BETA   =   DNRM2 ( M, U, 1 )

      IF (BETA .GT. ZERO) THEN
         CALL DSCAL ( M, (ONE / BETA), U, 1 )
         CALL APROD ( 2, M, N, V, U, LENIW, LENRW, IW, RW )
         ALFA   =   DNRM2 ( N, V, 1 )
      END IF

      IF (ALFA .GT. ZERO) THEN
         CALL DSCAL ( N, (ONE / ALFA), V, 1 )
         CALL DCOPY ( N, V, 1, W, 1 )
      END IF

      ARNORM =   ALFA * BETA
      IF (ARNORM .EQ. ZERO) GO TO 800

      RHOBAR =   ALFA
      PHIBAR =   BETA
      BNORM  =   BETA
      RNORM  =   BETA

 90   IF (NOUT   .GT.  0  ) THEN
         IF (DAMPSQ .EQ. ZERO) THEN
             WRITE(NOUT, 1200)
         ELSE
             WRITE(NOUT, 1300)
         END IF
         TEST1  = ONE
         TEST2  = ALFA / BETA
         WRITE(NOUT, 1500) ITN, X(1), RNORM, TEST1, TEST2
         WRITE(NOUT, 1600)
      END IF

*     ------------------------------------------------------------------
*     Main iteration loop.
*     ------------------------------------------------------------------
  100 ITN    = ITN + 1
c slmod begin
c ----------------------------------------------------------------------
c...  Store previous x values:
      OLDX(1:N) = X(1:N)       
c ----------------------------------------------------------------------
c slmod end
*     Perform the next step of the bidiagonalization to obtain the
*     next  BETA, U, ALFA, V.  These satisfy the relations
*                BETA*U  =  A*V  -  ALFA*U,
*                ALFA*V  =  A(transpose)*U  -  BETA*V.

      CALL DSCAL ( M, (- ALFA), U, 1 )
      CALL APROD ( 1, M, N, V, U, LENIW, LENRW, IW, RW )
      BETA   =   DNRM2 ( M, U, 1 )
      BBNORM =   BBNORM  +  ALFA**2  +  BETA**2  +  DAMPSQ

      IF (BETA .GT. ZERO) THEN
         CALL DSCAL ( M, (ONE / BETA), U, 1 )
         CALL DSCAL ( N, (- BETA), V, 1 )
         CALL APROD ( 2, M, N, V, U, LENIW, LENRW, IW, RW )
         ALFA   =   DNRM2 ( N, V, 1 )
         IF (ALFA .GT. ZERO) THEN
            CALL DSCAL ( N, (ONE / ALFA), V, 1 )
         END IF
      END IF

*     Use a plane rotation to eliminate the damping parameter.
*     This alters the diagonal (RHOBAR) of the lower-bidiagonal matrix.

      RHBAR2 = RHOBAR**2  +  DAMPSQ
      RHBAR1 = SQRT( RHBAR2 )
      CS1    = RHOBAR / RHBAR1
      SN1    = DAMP   / RHBAR1
      PSI    = SN1 * PHIBAR
      PHIBAR = CS1 * PHIBAR

*     Use a plane rotation to eliminate the subdiagonal element (BETA)
*     of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.

      RHO    =   SQRT( RHBAR2  +  BETA**2 )
      CS     =   RHBAR1 / RHO
      SN     =   BETA   / RHO
      THETA  =   SN * ALFA
      RHOBAR = - CS * ALFA
      PHI    =   CS * PHIBAR
      PHIBAR =   SN * PHIBAR
      TAU    =   SN * PHI

*     Update  X, W  and the standard error estimates.

      T1     =   PHI   / RHO
      T2     = - THETA / RHO
      T3     =   ONE   / RHO

      DO 200  I =  1, N
         T      =  W(I)
         X(I)   =  T1*T  +  X(I)
         W(I)   =  T2*T  +  V(I)
         T      = (T3*T)**2
         SE(I)  =  T     +  SE(I)
         DDNORM =  T     +  DDNORM
  200 CONTINUE

*     Use a plane rotation on the right to eliminate the
*     super-diagonal element (THETA) of the upper-bidiagonal matrix.
*     Then use the result to estimate  norm(X).

      DELTA  =   SN2 * RHO
      GAMBAR = - CS2 * RHO
      RHS    =   PHI    - DELTA * Z
      ZBAR   =   RHS    / GAMBAR
      XNORM  =   SQRT( XXNORM    + ZBAR **2 )
      GAMMA  =   SQRT( GAMBAR**2 + THETA**2 )
      CS2    =   GAMBAR / GAMMA
      SN2    =   THETA  / GAMMA
      Z      =   RHS    / GAMMA
      XXNORM =   XXNORM + Z**2

*     Test for convergence.
*     First, estimate the norm and condition of the matrix  Abar,
*     and the norms of  rbar  and  Abar(transpose)*rbar.

      ANORM  =   SQRT( BBNORM )
      ACOND  =   ANORM * SQRT( DDNORM )
      RES1   =   PHIBAR**2
      RES2   =   RES2  +  PSI**2
      RNORM  =   SQRT( RES1 + RES2 )
      ARNORM =   ALFA  * ABS( TAU )

*     Now use these norms to estimate certain other quantities,
*     some of which will be small near a solution.

      TEST1  =   RNORM /  BNORM
      TEST2  =   ZERO
      IF (RNORM .GT. ZERO) TEST2 = ARNORM / (ANORM * RNORM)
      TEST3  =   ONE   /  ACOND
      T1     =   TEST1 / (ONE  +  ANORM * XNORM / BNORM)
      RTOL   =   BTOL  +  ATOL *  ANORM * XNORM / BNORM

*     The following tests guard against extremely small values of
*     ATOL, BTOL  or  CTOL.  (The user may have set any or all of
*     the parameters  ATOL, BTOL, CONLIM  to zero.)
*     The effect is equivalent to the normal tests using
*     ATOL = RELPR,  BTOL = RELPR,  CONLIM = 1/RELPR.

      T3     =   ONE + TEST3
      T2     =   ONE + TEST2
      T1     =   ONE + T1
      IF (ITN .GE. ITNLIM) ISTOP = 7
      IF (T3  .LE. ONE   ) ISTOP = 6
      IF (T2  .LE. ONE   ) ISTOP = 5
      IF (T1  .LE. ONE   ) ISTOP = 4

*     Allow for tolerances set by the user.

      IF (TEST3 .LE. CTOL) ISTOP = 3
      IF (TEST2 .LE. ATOL) ISTOP = 2
      IF (TEST1 .LE. RTOL) ISTOP = 1
c slmod begin
c ----------------------------------------------------------------------
c...  Check if solution has a negative value:
      IF (STOP8OKAY) THEN
        DO I = 1, N
          IF (X(I).LT.0.0D0) THEN
            ISTOP = 8
            EXIT
          ENDIF
        ENDDO
      ENDIF
c ----------------------------------------------------------------------
c slmod end
*     ==================================================================

*     See if it is time to print something.

      IF (NOUT  .LE.  0       ) GO TO 600
      IF (N     .LE. 40       ) GO TO 400
      IF (ITN   .LE. 10       ) GO TO 400
      IF (ITN   .GE. ITNLIM-10) GO TO 400
      IF (MOD(ITN,10) .EQ. 0  ) GO TO 400
      IF (TEST3 .LE.  2.0*CTOL) GO TO 400
      IF (TEST2 .LE. 10.0*ATOL) GO TO 400
      IF (TEST1 .LE. 10.0*RTOL) GO TO 400
      IF (ISTOP .NE.  0       ) GO TO 400
      GO TO 600

*     Print a line for this iteration.

  400 WRITE(NOUT, 1500) ITN, X(1), RNORM, TEST1, TEST2, ANORM, ACOND
      IF (MOD(ITN,10) .EQ. 0) WRITE(NOUT, 1600)
*     ==================================================================

*     Stop if appropriate.
*     The convergence criteria are required to be met on  NCONV
*     consecutive iterations, where  NCONV  is set below.
*     Suggested value:  NCONV = 1, 2  or  3.

  600 IF (ISTOP .EQ. 0) NSTOP = 0
      IF (ISTOP .EQ. 0) GO TO 100
      NCONV  =   1
      NSTOP  =   NSTOP + 1
      IF (NSTOP .LT. NCONV  .AND.  ITN .LT. ITNLIM) ISTOP = 0
      IF (ISTOP .EQ. 0) GO TO 100
*     ------------------------------------------------------------------
*     End of iteration loop.
*     ------------------------------------------------------------------
c slmod begin
c ----------------------------------------------------------------------
c...  Restore last non-negative x solution:
      IF (ISTOP.EQ.8) X(1:N) = OLDX(1:N)
c ----------------------------------------------------------------------
c slmod end
*     Finish off the standard error estimates.

      T    =   ONE
      IF (M      .GT.   N )  T = M - N
      IF (DAMPSQ .GT. ZERO)  T = M
      T    =   RNORM / SQRT( T )

      DO 700  I = 1, N
         SE(I)  = T * SQRT( SE(I) )
  700 CONTINUE

*     Print the stopping condition.

  800 IF (NOUT .GT. 0) THEN
         WRITE(NOUT, 2000) EXIT, ISTOP, ITN,
     $                     EXIT, ANORM, ACOND,
     $                     EXIT, RNORM, ARNORM,
     $                     EXIT, BNORM, XNORM
         WRITE(NOUT, 3000) EXIT, MSG(ISTOP)
      END IF

  900 RETURN

*     ------------------------------------------------------------------
 1000 FORMAT(// 1P, A, '  Least-squares solution of  A*x = b'
     $    / ' The matrix  A  has', I7, ' rows   and', I7, ' columns'
     $    / ' The damping parameter is         DAMP   =', E10.2
     $    / ' ATOL   =', E10.2, 15X,        'CONLIM =', E10.2
     $    / ' BTOL   =', E10.2, 15X,        'ITNLIM =', I10)
 1200 FORMAT(// '   Itn       x(1)           Function',
     $   '     Compatible   LS        Norm A    Cond A' /)
 1300 FORMAT(// '   Itn       x(1)           Function',
     $   '     Compatible   LS     Norm Abar Cond Abar' /)
 1500 FORMAT(1P, I6, 2E17.9, 4E10.2)
 1600 FORMAT(1X)
 2000 FORMAT(/ 1P, A, 6X, 'ISTOP =', I3,   16X, 'ITN    =', I9
     $       /     A, 6X, 'ANORM =', E13.5, 6X, 'ACOND  =', E13.5
     $       /     A, 6X, 'RNORM =', E13.5, 6X, 'ARNORM =', E13.5,
     $       /     A, 6X, 'BNORM =', E13.5, 6X, 'XNORM  =', E13.5)
 3000 FORMAT( A, 6X, A )
*     ------------------------------------------------------------------
*     End of LSQR
      END
C     ******************************************************
C
C     WARNING.  Delete the following imitation BLAS routines
C               if a genuine BLAS library is available.
C
C     ******************************************************
C
      SUBROUTINE DCOPY ( N,X,INCX,Y,INCY )
      implicit none
! jdemod - add variables that weren't explicitly declared
      integer :: i
      INTEGER            N,INCX,INCY
      DOUBLE PRECISION   X(N),Y(N)
C
C     This may be replaced by the corresponding  BLAS  routine.
C     The following is a simple version for use with  LSQR.
C
      DO 10 I = 1, N
         Y(I) = X(I)
   10 CONTINUE
      RETURN
C
C     END OF DCOPY
      END
      DOUBLE PRECISION   FUNCTION DNRM2 ( N,X,INCX )
      INTEGER            N,INCX
      DOUBLE PRECISION   X(N)
C
C     This may be replaced by the corresponding  BLAS  routine.
C     The following is a simple version for use with  LSQR.
C
      INTEGER            I
      DOUBLE PRECISION   D, DSQRT
C
      D     = 0.0
      DO 10 I = 1, N
         D    = D + X(I)**2
   10 CONTINUE
      DNRM2 = DSQRT(D)
      RETURN
C
C     END OF DNRM2
      END
      SUBROUTINE DSCAL ( N,A,X,INCX )
      implicit none
! jdemod - add variables that weren't explicitly declared
      integer :: i
      INTEGER            N,INCX
      DOUBLE PRECISION   A,X(N)
C
C     This may be replaced by the corresponding  BLAS  routine.
C     The following is a simple version for use with  LSQR.
C
      DO 10 I = 1, N
         X(I) = A*X(I)
   10 CONTINUE
      RETURN
C
C     END OF DSCAL
      END
*********************************************************
*
*     These routines are for testing  LSQR.
*
*********************************************************

      SUBROUTINE APROD ( MODE, M, N, X, Y, LENIW, LENRW, IW, RW )
c slmod begin
      IMPLICIT none
c slmod end
      INTEGER            MODE, M, N, LENIW, LENRW
      INTEGER            IW(LENIW)
      DOUBLE PRECISION   X(N), Y(M), RW(LENRW)

*     ------------------------------------------------------------------
*     This is the matrix-vector product routine required by  LSQR
*     for a test matrix of the form  A = HY*D*HZ.  The quantities
*     defining D, HY, HZ are in the work array RW, followed by a
*     work array W.  These are passed to APROD1 and APROD2 in order to
*     make the code readable.
*     ------------------------------------------------------------------

      INTEGER            LOCD, LOCHY, LOCHZ, LOCW
c slmod begin
c ----------------------------------------------------------------------
      INTEGER i,j,ix,iy,i1,i2,codec,mfile,nfile,fp,nmap,irow,ncol,
     .        idum1,idum2
      LOGICAL ascii
      REAL    version
      REAL*8  VAL,ddum1
      INTEGER, ALLOCATABLE :: rind(:),cind(:)
      REAL*8 , ALLOCATABLE :: A(:,:), aval(:)
      CHARACTER file*(*),dummy*1024,cdum1*32

c      REAL*8  VAL,A(1000,1000)
c      LOGICAL FIRST
c      DATA FIRST /.TRUE./
      SAVE


c      IF (FIRST) THEN
c        DO I = 1, M
c          A(I,I) = DBLE(I)
c        ENDDO
c        WRITE(0,*) 'A:'
c        DO I = 1, M
c          WRITE(0,'(20F10.5)') A(I,1:N)
c        ENDDO
c        FIRST = .FALSE.
c      ENDIF

*                CALL APROD ( mode, m, n, x, y, LENIW, LENRW, IW, RW )
*
*     which must perform the following functions:
*
*                If MODE = 1, compute  y = y + A*x.
*                If MODE = 2, compute  x = x + A(transpose)*y.

c could use some vector power here...

      IF (codec.EQ.1) THEN
c...    'A' matix in full:
        SELECT CASE (mode)
          CASE (1)
            DO i = 1, m
              val = 0.0D0
              DO j = 1, n
                val = val + A(i,j) * x(j)
              ENDDO
              y(i) = y(i) + val
            ENDDO           
          CASE (2)
            DO i = 1, n
              val = 0.0D0
              DO j = 1, m
                val = val + A(j,i) * y(j)
              ENDDO
              x(i) = x(i) + val
            ENDDO           
        END SELECT

      ELSEIF (codec.EQ.2) THEN
c...    'A' matix is large and (hopefully) sparse:
        SELECT CASE (mode)
          CASE (1)
c            DO i = 1, m
c              val = 0.0D0
c              IF (rnum(i).GT.0) THEN
c                DO i1 = aind(i), aind(i)+rnum(i)-1
c                  j = cind(i1)
c                  val = val + aval(i1) * x(j)
c                ENDDO
c              ENDIF
c              y(i) = y(i) + val
c            ENDDO           
            DO i1 = 1, nmap
              i = rind(i1)
              j = cind(i1)
              y(i) = y(i) + aval(i1) * x(j)
            ENDDO
          CASE (2)
            DO i1 = 1, nmap
              j = rind(i1)
              i = cind(i1)
              x(i) = x(i) + aval(i1) * y(j)
            ENDDO
        END SELECT

      ENDIF

      RETURN



      ENTRY APROD_LOAD(m,n,file)
c      ENTRY APROD_LOAD(m,n,AMATRIX,load_mode)


c...  The A matrix should be loaded here, not in Main989 (although 
c     there should be the option to do that), since it makes no sense to be passing
c     around a potentially huge matrix, especially since passing it around might
c     be impossible once things get serious... LEFT OFF      

      fp = 99

      ascii = .TRUE.
      OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .     FORM='FORMATTED',STATUS='OLD',ERR=98)
      READ(fp,'(A1024)') dummy
      WRITE(0,*) 'dummy:',dummy(1:10)
      IF (dummy(1:1).EQ.'*') THEN
        BACKSPACE fp
      ELSE
c...    Try binary format:
        CLOSE(fp)
        ascii = .FALSE.
        OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .       FORM='UNFORMATTED',STATUS='OLD',ERR=98)   
        READ(fp) idum1
        IF (idum1.NE.888888888) 
     .    CALL ER('Aprod_Load','Unable to open inversion map file',*99)
      ENDIF


c      OPEN(fp,FILE=file(1:LEN_TRIM(file)),             ! *** Put this in a subroutine ***
c     .     FORM='FORMATTED',STATUS='OLD',ERR=98)
      IF (ascii) THEN
        dummy(1:1) = '*'
        DO WHILE (dummy(1:1).EQ.'*')
          READ(fp,'(A1024)') dummy
        ENDDO
        BACKSPACE fp
        READ(fp,'(A14,F5.1)') dummy,version
        IF (version.EQ.1.0) THEN
          READ(fp,'(A14,I7)') dummy,codec
          READ(fp,'(A14,I7)') dummy,mfile
          READ(fp,'(A14,I7)') dummy,nfile
          READ(fp,*)
          READ(fp,*)
          READ(fp,*)
        ELSE
          CALL ER('Aprod_Load','Unrecognized ASCII .map version',*99)
        ENDIF
      ELSE
        READ(fp) version
        IF (version.EQ.1.0) THEN
          READ(fp) codec
          READ(fp) mfile,nfile,idum1,(idum2,i1=1,idum1*2)
        ELSE
          CALL ER('Aprod_Load','Unrecognized binary .map version',*99)
        ENDIF
      ENDIF

      IF (mfile.NE.m.OR.nfile.NE.n) 
     .  CALL ER('AProd_Load','File A matix rank does not match '//
     .          'rank in call to AProd_Load',*99)

      WRITE(0,*) 'CODEC=',codec
      WRITE(0,*) 'ASCII=',ascii

      IF     (codec.EQ.1) THEN

        ALLOCATE(A(m,n))
        DO ix = 1, m
          DO i1 = 1, n, 5  ! Better way to write this?
            READ(fp,'(10D20.12)') (A(ix,iy),iy=i1,MIN(i1+4,n))
          ENDDO
c          IF (n.LE.5) 
c            WRITE(0,'(A,I4,10D20.12)') 'A:',ix,(A(ix,iy),iy=1,n)
        ENDDO        

        IF (n.LE.6) THEN
          WRITE(0,*) 'A:'
          DO I = 1, m
            WRITE(0,'(10F10.5)') (a(i,j),j=1,n)
          ENDDO
        ENDIF

      ELSEIF (codec.EQ.2) THEN
c...    Scan file to determine the array sizes required:

        WRITE(0,*) 'Scanning .map file...'
        nmap = 0
        DO WHILE (.TRUE.) 
          IF (ascii) THEN
            READ(fp,'(A1024)',END=10) dummy
            IF (dummy(1:5).EQ.'pixel') THEN
              READ(dummy,'(A6,2I8)') cdum1,irow,ncol
              nmap = nmap + ncol
            ENDIF
          ELSE
            READ(fp,END=10) idum1
            IF     (idum1.EQ.777777777) THEN
              READ(fp) irow,ncol,(idum1,ddum1,i1=1,ncol)
              nmap = nmap + ncol
            ELSEIF (idum1.EQ.999999999) THEN
              EXIT
            ENDIF
          ENDIF
        ENDDO
 10     CONTINUE
        WRITE(0,*) 'NMAP:',nmap
        WRITE(0,'(1X,A,F10.1,A)') 
     .    'ASKING FOR ',REAL((4+4+8)*nmap)/1.E6,' MBYTES'

        ALLOCATE(rind(nmap))
        ALLOCATE(cind(nmap))
        ALLOCATE(aval(nmap))

        nmap = 0
        REWIND(fp)
        DO WHILE (.TRUE.) 
          IF (ascii) THEN
            READ(fp,'(A1024)',END=20) dummy
            IF (dummy(1:5).EQ.'pixel') THEN
              READ(dummy,'(A6,2I8)') cdum1,irow,ncol

c              IF (nmap+ncol.GT.2115748)THEN
c              IF (irow.EQ.528191)THEN
c                WRITE(0,*) 'DUMMY: '//dummy(1:50)
c                WRITE(0,*) 'DUMMY: ',nmap,ncol
c                STOP 'GOTCHA!'
c              ENDIF

              rind(nmap+1:nmap+ncol) = irow 
              READ(fp,'(4(I7,D22.15))') 
     .          (cind(nmap+i1),aval(nmap+i1),i1=1,ncol)
              nmap = nmap + ncol             
            ENDIF
          ELSE
            READ(fp,END=20) idum1
            IF (idum1.EQ.777777777) THEN
              READ(fp) irow,ncol,(cind(nmap+i1),aval(nmap+i1),i1=1,ncol)
              rind(nmap+1:nmap+ncol) = irow 
              nmap = nmap + ncol             
            ENDIF
          ENDIF
        ENDDO
 20     CONTINUE
        WRITE(0,*) 'NMAP:',nmap

        IF (n.LE.6) THEN
          WRITE(0,*) 'A:'
          DO i1 = 1, nmap
            WRITE(0,'(2I7,F10.5)') rind(i1),cind(i1),aval(i1)
          ENDDO
        ENDIF

c        STOP 'sdfsd'

      ELSE
        CALL ER('AProd_Load','Bad codec in .map file',*99)
      ENDIF

      CLOSE(fp)



      RETURN


      ENTRY APROD_CLOSE

      IF     (codec.EQ.1) THEN
        WRITE(0,*) '  DEALLOCATING A'
        IF (ALLOCATED(A)) DEALLOCATE(A)
        WRITE(0,*) '  DONE'
      ELSEIF (codec.EQ.2) THEN
        WRITE(0,*) '  DEALLOCATING RIND,CIND,AVAL'
        IF (ALLOCATED(rind)) DEALLOCATE(rind)
c        WRITE(0,*) '  DEALLOCATING RIND'
        IF (ALLOCATED(cind)) DEALLOCATE(cind)
c        WRITE(0,*) '  DEALLOCATING AVAL'
        IF (ALLOCATED(aval)) DEALLOCATE(aval)
        WRITE(0,*) '  DONE'
      ELSE
        CALL ER('AProd','Unknown LOAD_MODE for close',*99)
      ENDIF

      RETURN
 98   CALL ER('Main989','Unable to open inversion map file',*99)
 99   STOP
c ----------------------------------------------------------------------
c slmod end
      LOCD   = 1
      LOCHY  = LOCD  + N
      LOCHZ  = LOCHY + M
      LOCW   = LOCHZ + N

      IF (MODE .EQ. 1) CALL APROD1( M, N, X, Y,
     $   RW(LOCD), RW(LOCHY), RW(LOCHZ), RW(LOCW) )

      IF (MODE .NE. 1) CALL APROD2( M, N, X, Y,
     $   RW(LOCD), RW(LOCHY), RW(LOCHZ), RW(LOCW) )

*     End of APROD
      END

      SUBROUTINE APROD1( M, N, X, Y, D, HY, HZ, W )
      implicit none
      INTEGER            M, N
      DOUBLE PRECISION   X(N), Y(M), D(N), HY(M), HZ(N), W(M)
      
*     ------------------------------------------------------------------
*     APROD1  computes  Y = Y + A*X  for subroutine APROD,
*     where A is a test matrix of the form  A = HY*D*HZ,
*     and the latter matrices HY, D, HZ are represented by
*     input vectors with the same name.
*     ------------------------------------------------------------------

      INTEGER            I
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0 )

      CALL HPROD ( N, HZ, X, W )

      DO 100 I = 1, N
         W(I)  = D(I) * W(I)
  100 CONTINUE

      DO 200 I = N + 1, M
         W(I)  = ZERO
  200 CONTINUE
      
      CALL HPROD ( M, HY, W, W )

      DO 600 I = 1, M
         Y(I)  = Y(I) + W(I)
  600 CONTINUE

*     End of APROD1
      END

      SUBROUTINE APROD2( M, N, X, Y, D, HY, HZ, W )
      implicit none
      INTEGER            M, N
      DOUBLE PRECISION   X(N), Y(M), D(N), HY(M), HZ(N), W(M)

*     ------------------------------------------------------------------
*     APROD2  computes  X = X + A(T)*Y  for subroutine APROD,
*     where  A  is a test matrix of the form  A = HY*D*HZ,
*     and the latter matrices  HY, D, HZ  are represented by
*     input vectors with the same name.
*     ------------------------------------------------------------------

      INTEGER            I

      CALL HPROD ( M, HY, Y, W )

      DO 100 I = 1, N
         W(I)  = D(I)*W(I)
  100 CONTINUE

      CALL HPROD ( N, HZ, W, W )

      DO 600 I = 1, N
         X(I)  = X(I) + W(I)
  600 CONTINUE

*     End of APROD2
      END

      SUBROUTINE HPROD ( N, HZ, X, Y )
      implicit none
      INTEGER            N
      DOUBLE PRECISION   HZ(N), X(N), Y(N)

*     ------------------------------------------------------------------
*     HPROD  applies a Householder transformation stored in  HZ
*     to get  Y = ( I - 2*HZ*HZ(transpose) ) * X.
*     ------------------------------------------------------------------

      INTEGER            I
      DOUBLE PRECISION   S

      S      = 0.0
      DO 100 I = 1, N
         S     = HZ(I) * X(I)  +  S
  100 CONTINUE

      S      = S + S
      DO 200 I = 1, N
         Y(I)  = X(I)  -  S * HZ(I)
  200 CONTINUE

*     End of HPROD
      END

      SUBROUTINE LSTP  ( M, N, NDUPLC, NPOWER, DAMP, X,
     $                   B, D, HY, HZ, W, ACOND, RNORM )
      implicit none
      INTEGER            M, N, MAXMN, NDUPLC, NPOWER
      DOUBLE PRECISION   DAMP, ACOND, RNORM
      DOUBLE PRECISION   B(M), X(N), D(N), HY(M), HZ(N), W(M)

*     ------------------------------------------------------------------
*     LSTP  generates a sparse least-squares test problem of the form
*
*                (   A    )*X = ( B ) 
*                ( DAMP*I )     ( 0 )
*
*     having a specified solution X.  The matrix A is constructed
*     in the form  A = HY*D*HZ,  where D is an M by N diagonal matrix,
*     and HY and HZ are Householder transformations.
*
*     The first 6 parameters are input to LSTP.  The remainder are
*     output.  LSTP is intended for use with M .GE. N.
*
*
*     Functions and subroutines
*
*     TESTPROB           APROD1, HPROD
*     BLAS               DNRM2
*     ------------------------------------------------------------------

*     Intrinsics and local variables

      INTRINSIC          COS,  SIN, SQRT
      INTEGER            I, J
      DOUBLE PRECISION   DNRM2
      DOUBLE PRECISION   ALFA, BETA, DAMPSQ, FOURPI, T
      DOUBLE PRECISION   ZERO,        ONE
      PARAMETER        ( ZERO = 0.0,  ONE = 1.0 )

*     ------------------------------------------------------------------
*     Make two vectors of norm 1.0 for the Householder transformations.
*     FOURPI  need not be exact.
*     ------------------------------------------------------------------
      DAMPSQ = DAMP**2
      FOURPI = 4.0 * 3.141592
      ALFA   = FOURPI / M
      BETA   = FOURPI / N

      DO 100 I = 1, M
         HY(I) = SIN( I * ALFA )
  100 CONTINUE

      DO 200 I = 1, N
         HZ(I) = COS( I * BETA )
  200 CONTINUE                

      ALFA   = DNRM2 ( M, HY, 1 )
      BETA   = DNRM2 ( N, HZ, 1 )
      CALL DSCAL ( M, (- ONE / ALFA), HY, 1 )
      CALL DSCAL ( N, (- ONE / BETA), HZ, 1 )
*            
*     ------------------------------------------------------------------
*     Set the diagonal matrix  D.  These are the singular values of  A.
*     ------------------------------------------------------------------
      DO 300 I = 1, N
         J     = (I - 1 + NDUPLC) / NDUPLC
         T     =  J * NDUPLC
         T     =  T / N
         D(I)  =  T**NPOWER
  300 CONTINUE

      ACOND  = SQRT( (D(N)**2 + DAMPSQ) / (D(1)**2 + DAMPSQ) )

*     ------------------------------------------------------------------
*     Compute the residual vector, storing it in  B.
*     It takes the form  HY*( s )
*                           ( t )
*     where  s  is obtained from  D*s = DAMP**2 * HZ * X
*     and    t  can be anything.
*     ------------------------------------------------------------------
      CALL HPROD ( N, HZ, X, B )

      DO 500 I = 1, N
         B(I)  = DAMPSQ * B(I) / D(I)
  500 CONTINUE

      T      = ONE
      DO 600 I =   N + 1, M
         J     =   I - N
         B(I)  =  (T * J) / M
         T     = - T
  600 CONTINUE

      CALL HPROD ( M, HY, B, B )

*     ------------------------------------------------------------------
*     Now compute the true  B  =  RESIDUAL  +  A*X.
*     ------------------------------------------------------------------
      RNORM  = SQRT(            DNRM2 ( M, B, 1 )**2
     $              +  DAMPSQ * DNRM2 ( N, X, 1 )**2 )
      CALL APROD1( M, N, X, B, D, HY, HZ, W )

*     End of LSTP
      END

      SUBROUTINE TEST_LSQR  ( M, N, NDUPLC, NPOWER, DAMP )
      implicit none
! jdemod - add variables that weren't explicitly declared
!          they are parameters
      integer :: maxm,maxn,leniw,lenrw, i, itn

      INTEGER            M, N, NDUPLC, NPOWER
      DOUBLE PRECISION   DAMP

*     ------------------------------------------------------------------
*     This is an example driver routine for running LSQR.
*     It generates a test problem, solves it, and examines the results.
*     Note that subroutine APROD must be declared EXTERNAL
*     if it is used only in the call to LSQR.
*
*
*     Functions and subroutines
*
*     TESTPROB           APROD
*     BLAS               DCOPY, DNRM2, DSCAL
*     ------------------------------------------------------------------

*     Intrinsics and local variables

      INTRINSIC          MAX, SQRT
      EXTERNAL           APROD
      INTEGER            ISTOP, ITNLIM, J, NOUT
      DOUBLE PRECISION   DNRM2
    
      PARAMETER        ( MAXM = 200,  MAXN = 100 )
      DOUBLE PRECISION   B(MAXM),  U(MAXM),
     $                   V(MAXN),  W(MAXN), X(MAXN),
     $                   SE(MAXN), XTRUE(MAXN)
      DOUBLE PRECISION   ATOL, BTOL, CONLIM,
     $                   ANORM, ACOND, RNORM, ARNORM,
     $                   DAMPSQ, ENORM, ETOL, XNORM

      PARAMETER        ( LENIW = 1,  LENRW = 600 )
      INTEGER            IW(LENIW)
      DOUBLE PRECISION   RW(LENRW)
      INTEGER            LOCD, LOCHY, LOCHZ, LOCW, LTOTAL

      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE = 1.0 )

      CHARACTER*34       LINE
      DATA               LINE
     $                 /'----------------------------------'/
c slmod begin
c ----------------------------------------------------------------------

      XTRUE = 0.0D0
      XTRUE(1:M) = 1.0D0

      DO I = 1, M
        U(I) = 1.0/DBLE(I)
      ENDDO

      WRITE(0,*) 'b:'
      DO I = 1, M
        WRITE(0,'(F10.5)') U(I)
      ENDDO

      DAMP = 0.0D0

      ATOL   = 1.0D-10
      BTOL   = ATOL
      CONLIM = 10.0 * ACOND
      ITNLIM = M + N + 50
      NOUT   = 6



      CALL LSQR  ( M, N, APROD, DAMP,
     $             LENIW, LENRW, IW, RW,
     $             U, V, W, X, SE,
     $             ATOL, BTOL, CONLIM, ITNLIM, NOUT,
     $             ISTOP, ITN, ANORM, ACOND, RNORM, ARNORM, XNORM )




      WRITE(0,*) 'x:',istop
      DO I = 1, MIN(20,N)
        WRITE(0,'(F10.5)') X(I)
      ENDDO

      DAMPSQ = DAMP**2
      WRITE(NOUT, 2000)
      WRITE(NOUT, 2100) RNORM, ARNORM, XNORM

*     Compute  U = A*X - B.
*     This is the negative of the usual residual vector.
*     It will be close to zero only if  B  is a compatible rhs
*     and  X  is close to a solution.

      CALL DCOPY ( M, B, 1, U, 1 )
      CALL DSCAL ( M, (-ONE), U, 1 )
      CALL APROD ( 1, M, N, X, U, LENIW, LENRW, IW, RW )

      WRITE(0,*) 'u:',m
      DO I = 1, M
        WRITE(0,'(F10.5)') u(I)
      ENDDO

*     Compute  V = A(transpose)*U  +  DAMP**2 * X.
*     This will be close to zero in all cases
*     if  X  is close to a solution.

      CALL DCOPY ( N, X, 1, V, 1 )
      CALL DSCAL ( N, DAMPSQ, V, 1 )
      CALL APROD ( 2, M, N, V, U, LENIW, LENRW, IW, RW )

*     Compute the norms associated with  X, U, V.

      XNORM  = DNRM2 ( N, X, 1 )
      RNORM  = SQRT( DNRM2 ( M, U, 1 )**2  +  DAMPSQ * XNORM**2 )
      ARNORM = DNRM2 ( N, V, 1 )
      WRITE(NOUT, 2200) RNORM, ARNORM, XNORM

*     Print the solution and standard error estimates from  LSQR.

      WRITE(NOUT, 2500) (J, X(J),  J = 1, N)
      WRITE(NOUT, 2600) (J, SE(J), J = 1, N)

*     Print a clue about whether the solution looks OK.
                 
      DO J = 1, N
         W(J)  = X(J) - XTRUE(J)
      ENDDO
      ENORM    = DNRM2 ( N, W, 1 ) / (ONE  +  DNRM2 ( N, XTRUE, 1 ))
      ETOL     = 1.0D-5
      IF (ENORM .LE. ETOL) WRITE(NOUT, 3000) ENORM
      IF (ENORM .GT. ETOL) WRITE(NOUT, 3100) ENORM



      RETURN

c ----------------------------------------------------------------------
c slmod end

*     Set the desired solution  XTRUE.

      DO 100 J = 1, N
         XTRUE(J) = 1.0
  100 CONTINUE

*     Generate the specified test problem.
*     The workspace array  IW  is not needed in this application.
*     The workspace array  RW  is used for the following vectors:
*     D(N), HY(M), HZ(N), W(MAX(M,N)).
*     The vectors  D, HY, HZ  will define the test matrix A.
*     W is needed for workspace in APROD1 and APROD2.

      LOCD   = 1
      LOCHY  = LOCD  + N
      LOCHZ  = LOCHY + M
      LOCW   = LOCHZ + N
      LTOTAL = LOCW  + MAX(M,N) - 1
      IF (LTOTAL .GT. LENRW) GO TO 900

      CALL LSTP  ( M, N, NDUPLC, NPOWER, DAMP, XTRUE,
     $             B, RW(LOCD), RW(LOCHY), RW(LOCHZ), RW(LOCW),
     $             ACOND, RNORM )

*     Solve the problem defined by APROD, DAMP and B.
*     Copy the rhs vector B into U  (LSQR will overwrite U)
*     and set the other input parameters for LSQR.

      CALL DCOPY ( M, B, 1, U, 1 )
      ATOL   = 1.0E-10
      BTOL   = ATOL
      CONLIM = 10.0 * ACOND
      ITNLIM = M + N + 50
      NOUT   = 6
      WRITE(NOUT, 1000) LINE, LINE,
     $                  M, N, NDUPLC, NPOWER, DAMP, ACOND, RNORM,
     $                  LINE, LINE

      CALL LSQR  ( M, N, APROD, DAMP,
     $             LENIW, LENRW, IW, RW,
     $             U, V, W, X, SE,
     $             ATOL, BTOL, CONLIM, ITNLIM, NOUT,
     $             ISTOP, ITN, ANORM, ACOND, RNORM, ARNORM, XNORM )

*     Examine the results.
*     We print the residual norms  RNORM  and  ARNORM  given by LSQR,
*     and then compute their true values in terms of the solution  X
*     obtained by  LSQR.  At least one of them should be small.

      DAMPSQ = DAMP**2
      WRITE(NOUT, 2000)
      WRITE(NOUT, 2100) RNORM, ARNORM, XNORM

*     Compute  U = A*X - B.
*     This is the negative of the usual residual vector.
*     It will be close to zero only if  B  is a compatible rhs
*     and  X  is close to a solution.

      CALL DCOPY ( M, B, 1, U, 1 )
      CALL DSCAL ( M, (-ONE), U, 1 )
      CALL APROD ( 1, M, N, X, U, LENIW, LENRW, IW, RW )

*     Compute  V = A(transpose)*U  +  DAMP**2 * X.
*     This will be close to zero in all cases
*     if  X  is close to a solution.

      CALL DCOPY ( N, X, 1, V, 1 )
      CALL DSCAL ( N, DAMPSQ, V, 1 )
      CALL APROD ( 2, M, N, V, U, LENIW, LENRW, IW, RW )

*     Compute the norms associated with  X, U, V.

      XNORM  = DNRM2 ( N, X, 1 )
      RNORM  = SQRT( DNRM2 ( M, U, 1 )**2  +  DAMPSQ * XNORM**2 )
      ARNORM = DNRM2 ( N, V, 1 )
      WRITE(NOUT, 2200) RNORM, ARNORM, XNORM

*     Print the solution and standard error estimates from  LSQR.

      WRITE(NOUT, 2500) (J, X(J),  J = 1, N)
      WRITE(NOUT, 2600) (J, SE(J), J = 1, N)

*     Print a clue about whether the solution looks OK.
                 
      DO 500 J = 1, N
         W(J)  = X(J) - XTRUE(J)
  500 CONTINUE
      ENORM    = DNRM2 ( N, W, 1 ) / (ONE  +  DNRM2 ( N, XTRUE, 1 ))
      ETOL     = 1.0D-5
      IF (ENORM .LE. ETOL) WRITE(NOUT, 3000) ENORM
      IF (ENORM .GT. ETOL) WRITE(NOUT, 3100) ENORM
      RETURN

*     Not enough workspace.

  900 WRITE(NOUT, 9000) LTOTAL
      RETURN
                                 
 1000 FORMAT(1P
     $ // 1X, 2A
     $ /  ' Least-Squares Test Problem      P(', 4I5, E12.2, ' )'
     $ // ' Condition no. =', E12.4,  '     Residual function =', E17.9
     $ /  1X, 2A)
 2000 FORMAT(
     $ // 22X, ' Residual norm    Residual norm    Solution norm'
     $  / 22X, '(Abar X - bbar)   (Normal eqns)         (X)' /)
 2100 FORMAT(1P, ' Estimated by LSQR', 3E17.5)
 2200 FORMAT(1P, ' Computed from  X ', 3E17.5) 
 2500 FORMAT(//' Solution  X' / 4(I6, G14.6))
 2600 FORMAT(/ ' Standard errors  SE' / 4(I6, G14.6))
 3000 FORMAT(1P / ' LSQR  appears to be successful.',
     $        '     Relative error in  X  =', E10.2)
 3100 FORMAT(1P / ' LSQR  appears to have failed.  ',
     $        '     Relative error in  X  =', E10.2)
 9000 FORMAT(/ ' XXX  Insufficient workspace.',
     $        '  The length of  RW  should be at least', I6)
*     End of TEST_LSQR
      END

      SUBROUTINE LSQR_MAIN
      implicit none
*     -------------
*     Main program.
*     -------------
      DOUBLE PRECISION   DAMP1, DAMP2, DAMP3, DAMP4, ZERO
*
      ZERO   = 0.0
      DAMP1  = 0.1
      DAMP2  = 0.01
      DAMP3  = 0.001
      DAMP4  = 0.0001
c slmod begin
      CALL TEST_LSQR  (  5, 10, 1, 1, ZERO  )
c
c      CALL TEST  (  1,  1, 1, 1, ZERO  )
c      CALL TEST  (  1,  1, 1, 1, ZERO  )
c      CALL TEST  (  2,  1, 1, 1, ZERO  )
c      CALL TEST  ( 40, 40, 4, 4, ZERO  )
c      CALL TEST  ( 40, 40, 4, 4, DAMP2 )
c      CALL TEST  ( 80, 40, 4, 4, DAMP2 )
c slmod end
      STOP

*     End of main program for testing LSQR
      END
c slmod begin
c ----------------------------------------------------------------------
      SUBROUTINE INVERT_LSQR  ( m, n, afile, x, b, mode , damp)
c      SUBROUTINE INVERT_LSQR  ( m, n, A, x, b, mode )
      IMPLICIT none

      INTEGER m,n,mode
      LOGICAL cont,nonnegative,negative
      REAL*8 ::  x(n),b(m)      
      CHARACTER afile*(*)
c      REAL*8 ::  A(m,n),x(n),b(m)      

      REAL*8 u(m),v(n),w(m),se(n)
      REAL*8, ALLOCATABLE :: oldx(:)

c      REAL*8, ALLOCATABLE :: u(:),v(:),w(:),se(:),oldx(:)

*     ------------------------------------------------------------------
*     This is an example driver routine for running LSQR.
*     It generates a test problem, solves it, and examines the results.
*     Note that subroutine APROD must be declared EXTERNAL
*     if it is used only in the call to LSQR.
*
*
*     Functions and subroutines
*
*     TESTPROB           APROD
*     BLAS               DCOPY, DNRM2, DSCAL
*     ------------------------------------------------------------------

*     Intrinsics and local variables

      INTRINSIC          MAX, SQRT
      EXTERNAL           APROD
      INTEGER            ISTOP, ITNLIM, J, NOUT, ITN, I,itntot,itnmax,
     .                   i1,COUNT
      DOUBLE PRECISION   DNRM2
    
c      PARAMETER        ( MAXM = 200,  MAXN = 100 )
c      DOUBLE PRECISION   B(MAXM),  U(MAXM),
c     $                   V(MAXN),  W(MAXN), X(MAXN),
c     $                   SE(MAXN), XTRUE(MAXN)
      REAL*8             ATOL, BTOL, CONLIM, DAMP, 
     $                   ANORM, ACOND, RNORM, ARNORM,
     $                   DAMPSQ, ENORM, ETOL, XNORM,lsqr_r



      INTEGER    LENIW    , LENRW
      PARAMETER (LENIW = 1, LENRW = 600 )
c      INTEGER, ALLOCATABLE :: IW(:)
c      REAL*8 , ALLOCATABLE :: RW(:)
      INTEGER            IW(LENIW)
      DOUBLE PRECISION   RW(LENRW)

      INTEGER            LOCD, LOCHY, LOCHZ, LOCW, LTOTAL

c      DOUBLE PRECISION   ONE
c      PARAMETER        ( ONE = 1.0 )

c      CHARACTER*34       LINE
c      DATA               LINE
c     $                 /'----------------------------------'/

c      XTRUE = 0.0D0
c      XTRUE(1:M) = 1.0D0

      ALLOCATE(oldx(n))


c      ALLOCATE(u(m))
c      ALLOCATE(v(n))
c      ALLOCATE(w(m))
c      ALLOCATE(se(n))
c      ALLOCATE(iw(leniw))
c      ALLOCATE(rw(lenrw))

c      DO I = 1, M
c        U(I) = 1.0/DBLE(I)
c      ENDDO

      IF (mode.EQ.1) THEN
        nonnegative = .TRUE.
      ELSE
        nonnegative = .FALSE.
      ENDIF

      WRITE(0,*) 'NONNEGATIVE:',nonnegative

      u(1:m) = b(1:m)

      IF (m.LE.5.AND.n.LT.500) THEN
        WRITE(0,*) 'b:'
        DO I = 1, M
          WRITE(0,'(F10.5)') b(I)
        ENDDO
      ENDIF

      CALL APROD_Load(m,n,afile)

      itnmax = 100 * (M + N + 50)   ! An upper limit on the number of iterations

      itntot = 0
c...

      count = 0

      WRITE(0,*) 'LSQR DAMP:',damp

      IF (nonnegative) THEN

c        DAMP = 0.1D0
        ACOND = 1.0D0

        ATOL   = 1.0D-08          ! An estimate of the relative error in the data in A
        BTOL   = ATOL             ! An extimate of the relative error in the data in b
        CONLIM = 1.0D+09 * ACOND  ! An upper limit on cond(Abar)
        itnlim = -itnmax
        NOUT   = 6                ! File number for printed output

        u(1:m) = b(1:m)

        CALL LSQR  ( M, N, APROD, DAMP,
     $               LENIW, LENRW, IW, RW,
     $               U, V, W, X, SE,
     $               ATOL, BTOL, CONLIM, ITNLIM, NOUT,
     $               ISTOP, ITN, ANORM, ACOND, RNORM, ARNORM, XNORM )

        IF (m.LE.20.AND.n.LT.500) THEN
          WRITE(0,*) 'x:',istop,count,itn
          DO I = 1, N
            WRITE(0,'(2F14.9)') X(I),se(i)
          ENDDO
        ENDIF          

        itnlim = itnmax

c        DAMP = 0.1D0

        IF (istop.EQ.8) THEN


          DO WHILE(.TRUE.) 
 
            oldx(1:n) = x(1:n)
 
            u = 0.0D0
            CALL APROD ( 1, M, N, X, u, LENIW, LENRW, IW, RW )            

            IF (m.LE.25.AND.m.LT.500) THEN
              WRITE(0,*) '--------------------------'
              WRITE(0,*) 'Ax:'
              DO I = 1, m
                WRITE(0,'(2E20.10)') u(I)
              ENDDO
            ENDIF          

c...        Assign 'r':
            WRITE(0,*) '*** USING A + SIGN! (NOT WHAT REF. SAYS!) ***'

            u(1:m) = +(b(1:m) - u(1:m))
            
            IF (m.LE.25.AND.n.LT.500) THEN
              WRITE(0,*) '-r:'
              DO I = 1, m
                WRITE(0,'(2E20.10)') u(I)
              ENDDO
            ENDIF          

            lsqr_r = 0.0D0
            DO I = 1, m
              lsqr_r = lsqr_r + DSQRT(u(i)**2)
            ENDDO
            WRITE(0,*) '  LSQR_R:',lsqr_r

            CALL LSQR  ( M, N, APROD, DAMP,
     $                   LENIW, LENRW, IW, RW,
     $                   U, V, W, X, SE,
     $                   ATOL, BTOL, CONLIM, ITNLIM, NOUT,
     $                   ISTOP, ITN, ANORM, ACOND, RNORM, ARNORM, XNORM)

            IF (n.LE.25) THEN
              WRITE(0,*) 'w:',istop,count,itn
              DO I = 1, N
                WRITE(0,'(2E20.10)') X(I),se(i)
              ENDDO
            ENDIF          

c...        Update 'x':
            DO i1 = 1, n
              x(i1) = MAX(0.0D0,oldx(i1)+x(i1))
            ENDDO
           
            IF (n.LE.25) THEN
            WRITE(0,*) 'x:',istop,count,itn
              DO I = 1, N
                WRITE(0,'(2E20.10)') X(I)
              ENDDO
            ENDIF          

 
            count = count + 1
            WRITE(NOUT,*) 'NON-NEGATIVITY ITERATION:',count
            WRITE(0   ,*) 'NON-NEGATIVITY ITERATION:',count
            IF (count.EQ.10) EXIT

          ENDDO

        ENDIF

      ELSE


c        DAMP = 0.001D0 ! 0.01D0 ! 0.15D0 ! 0.1D0
        ACOND = 1.0D0

        ATOL   = 1.0D-08          ! An estimate of the relative error in the data in A
c        ATOL   = 1.0D-10          ! An estimate of the relative error in the data in A
        BTOL   = ATOL             ! An extimate of the relative error in the data in b
        CONLIM = 1.0D+09 * ACOND  ! An upper limit on cond(Abar)
        itnlim = itnmax
        NOUT   = 6                ! File number for printed output

        u(1:m) = b(1:m)

        CALL LSQR  ( M, N, APROD, DAMP,
     $               LENIW, LENRW, IW, RW,
     $               U, V, W, X, SE,
     $               ATOL, BTOL, CONLIM, ITNLIM, NOUT,
     $               ISTOP, ITN, ANORM, ACOND, RNORM, ARNORM, XNORM )

      ENDIF



      u = 0.0D0
      CALL APROD ( 1, M, N, X, u, LENIW, LENRW, IW, RW )

      
      IF (m.LE.20.AND.n.LT.500) THEN
        WRITE(0,*) 'x:'
        DO I = 1, N
          WRITE(0,'(2F14.9)') X(I),se(i)
        ENDDO
c        WRITE(0,*) 'b-Ax:'
c        DO I = 1, m
c          WRITE(0,'(F14.9)') b(i)-u(i)
c        ENDDO
      ELSE
c        WRITE(NOUT,*) 'x:'
c        DO I = 1, N
c          WRITE(NOUT,'(2E20.10)') X(I),se(i)
c        ENDDO
      ENDIF

c      IF (m.LE.5) THEN
c        WRITE(0,*) 'b:'
c        w = 0.0D0    
c        DO i = 1, m
c          DO j = 1, n
c            w(i) = w(i) + A(i,j) * x(j)
c          ENDDO
c          WRITE(0,'(F14.9)') w(I)
c        ENDDO
c      ENDIF

*     Compute  U = A*X - B.
*     This is the negative of the usual residual vector.
*     It will be close to zero only if  B  is a compatible rhs
*     and  X  is close to a solution.

      WRITE(0,*) 'CLOSING APROD'

      CALL APROD_Close

      WRITE(0,*) 'DEALLOCATING OLDX'

      IF (ALLOCATED(oldx)) DEALLOCATE(oldx)

c      WRITE(0,*) 'DEALLOCATING1'
c      IF (ALLOCATED(u)) DEALLOCATE(u)
c      WRITE(0,*) 'DEALLOCATING2'
c      IF (ALLOCATED(v)) DEALLOCATE(v)
c      WRITE(0,*) 'DEALLOCATING3'
c      IF (ALLOCATED(w)) DEALLOCATE(w)
c      WRITE(0,*) 'DEALLOCATING4'
c      IF (ALLOCATED(se)) DEALLOCATE(se)
c      WRITE(0,*) 'DEALLOCATING5'
c      DEALLOCATE(iw)
c      WRITE(0,*) 'DEALLOCATING6'
c      DEALLOCATE(rw)

c         STOP 'sdfsdfsd'


      WRITE(0,*) 'DONE IN INVERSION ROUTINE'

      RETURN
 99   STOP
      END
c
c ----------------------------------------------------------------------
c slmod end
