 
 
      SUBROUTINE EIRENE_LSBRR(NRA,NCA,A,LDA,B,TOL,X,RES,KBASIS)
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      INTEGER NRA,NCA,LDA,KBASIS,I,J,IAA
      REAL(DP) A(LDA,NCA),B(NRA),X(NCA),RES(NRA),TOL
      write (iunout,*) ' lsbrr is called EIRENE_'
      call EIRENE_exit_own(1)
      end
