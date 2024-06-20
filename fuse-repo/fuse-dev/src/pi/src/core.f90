!
! ======================================================================
!
SUBROUTINE piCore
#include "compile_opt.inc"
IMPLICIT none

INTEGER :: piCUDA,piOpenACC
INTEGER :: ierr


call piSerial

call piOpenMP

#if CUDA_FORTRAN == 1
ierr = piCuda()
#endif

ierr = piOpenACC()

RETURN  
END SUBROUTINE piCore
!
! ======================================================================
!
