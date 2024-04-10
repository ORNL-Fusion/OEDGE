!
! ======================================================================
!

! mpirun -np 4 ./pi

PROGRAM piMain
#include "compile_opt.inc"
USE mod_environment
USE mod_configure
USE mod_journal
IMPLICIT none

#if MPI == 1
WRITE(0,*) 'MPI activated'
#endif

CALL piInvestigateEnvironment

CALL piConfigureSimulation

CALL piInitializeJournal
  
! Build geometry

! Assign plasma
  
! Get rates

! Distribute setup data
  
! Generate random numbers

! Follow particles

CALL piCore
  
! Collect results

! Process results
  
END PROGRAM piMain
!
! ======================================================================
!
