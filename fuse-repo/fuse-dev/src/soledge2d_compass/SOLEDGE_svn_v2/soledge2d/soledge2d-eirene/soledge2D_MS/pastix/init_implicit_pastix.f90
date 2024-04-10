subroutine init_implicit_pastix(CSC_,nmat)
  use all_variables, only : global_parameters, zones
  use Mpastix_solve
  implicit none
#include "pastix_fortran.h"
  Type(CSC) :: CSC_
  integer*4 k,nmat
  real*8,allocatable,target :: dummy(:)
  integer :: OMP_GET_MAX_THREADS
  nmat=0
  do k=1,global_parameters%N_Zones
     nmat=nmat+(Zones(k)%mesh%Nx+2)*(Zones(k)%mesh%Nz+2)
  end do
  allocate(CSC_%colptr(1:nmat+1))
  allocate(CSC_%b(1:nmat))
  allocate(CSC_%loc2glb(1:nmat))
  allocate(CSC_%perm(1:nmat))
  allocate(CSC_%invp(1:nmat))

  allocate(dummy(1:1))
  CSC_%row => CSC_%colptr
  CSC_%avals => dummy

  CSC_%nmat=nmat
  CSC_%pastix_data = 0
  CSC_%pastix_comm = 0
  CSC_%rhs        = 1
  CSC_%iparm(IPARM_MODIFY_PARAMETER) = API_NO
  CSC_%iparm(IPARM_START_TASK)       = API_TASK_INIT
  CSC_%iparm(IPARM_END_TASK)         = API_TASK_INIT

  call pastix_fortran(CSC_%pastix_data,CSC_%pastix_comm, &
       CSC_%nmat, CSC_%colptr, CSC_%row, CSC_%avals, &
       CSC_%perm, CSC_%invp, CSC_%b, &
       CSC_%rhs, CSC_%iparm, CSC_%dparm)

  CSC_%iparm(IPARM_THREAD_NBR) = OMP_GET_MAX_THREADS()
  CSC_%iparm(IPARM_VERBOSE)          = API_VERBOSE_NOT !API_VERBOSE_YES
  CSC_%iparm(IPARM_SYM)              = API_SYM_NO
  CSC_%iparm(IPARM_FACTORIZATION)    = API_FACT_LU

  deallocate(dummy)
end subroutine init_implicit_pastix
