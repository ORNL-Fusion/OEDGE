

      SUBROUTINE LSBRR(NRA,NCA,A,LDA,B,TOL,X,RES,KBASIS)
      USE PRECISION
      IMPLICIT NONE
      INTEGER NRA,NCA,LDA,KBASIS,I,J,IAA
      REAL(DP) A(LDA,NCA),B(NRA),X(NCA),RES(NRA),TOL
      write (6,*) ' lsbrr is called '
      call exit_own(1)
      end
