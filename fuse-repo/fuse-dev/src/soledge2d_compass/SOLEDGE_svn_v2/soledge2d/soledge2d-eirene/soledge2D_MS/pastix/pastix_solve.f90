module Mpastix_solve

  implicit none

#include "pastix_fortran.h"

  Type :: CSC
     pastix_int_t,pointer :: colptr(:)
     pastix_int_t,pointer :: row(:)
     pastix_float_t,pointer :: avals(:)
     pastix_float_t,pointer :: b(:)
     pastix_int_t :: nmat,nnz
     pastix_data_ptr_t :: pastix_data
     pastix_data_ptr_t :: data_check
     integer :: pastix_comm
     pastix_int_t :: rhs
     pastix_int_t,allocatable :: loc2glb(:),perm(:),invp(:)
     pastix_int_t :: iparm(64)
     real*8 :: dparm(64)
  end type CSC

contains

  subroutine clear_pastix(CSC_)
    use Mlist
    use Msort
    implicit none
    Type(CSC) :: CSC_
    deallocate(CSC_%row)
    deallocate(CSC_%avals)
  end subroutine clear_pastix

  subroutine free_pastix(CSC_)
    use Mlist
    use Msort
    implicit none
    Type(CSC) :: CSC_
    deallocate(CSC_%colptr)
    deallocate(CSC_%b)
    deallocate(CSC_%loc2glb)
    deallocate(CSC_%perm)
    deallocate(CSC_%invp)
    CSC_%iparm(IPARM_START_TASK) = API_TASK_CLEAN
    CSC_%iparm(IPARM_END_TASK) = API_TASK_CLEAN
    call pastix_fortran(CSC_%pastix_data,CSC_%pastix_comm, &
         CSC_%nmat, CSC_%colptr, CSC_%row, CSC_%avals, &
         CSC_%perm, CSC_%invp, CSC_%b, &
         CSC_%rhs, CSC_%iparm, CSC_%dparm)
  end subroutine free_pastix

  subroutine fill_implicit_pastix(CSC_,nnz,matrix)
    use Mlist
    use Msort
    implicit none
    integer*4 nnz
    type(cell),intent(in), pointer :: matrix
    Type(CSC) :: CSC_
    integer*4 :: rows(nnz)
    integer*4 :: cols(nnz)
    integer*4 :: order(nnz)
    real*8 :: vals(nnz)
    integer*4 k,n
    pastix_int_t :: flagcor
    call list2array(rows,cols,vals,nnz,1,matrix)

!!$    open(unit=101,file='rows',status='unknown')
!!$    open(unit=102,file='cols',status='unknown')
!!$    open(unit=103,file='vals',status='unknown')
    do k=1,nnz
       order(k)=k
!!$       write(101,*) rows(k)
!!$       write(102,*) cols(k)
!!$       write(103,*) vals(k)
    end do
!!$    close(101)
!!$    close(102)
!!$    close(103)

    if(associated(CSC_%row)) then
       nullify(CSC_%row)
    end if
    if(associated(CSC_%avals)) then
       nullify(CSC_%avals)
    end if

    allocate(CSC_%row(1:nnz))
    allocate(CSC_%avals(1:nnz))

    call QsortC(cols,order)
    CSC_%nnz=nnz
    do k=1,nnz
       CSC_%row(k)=rows(order(k))
       CSC_%avals(k)=vals(order(k))
    end do
    CSC_%colptr(1)=1
    n=2
    do k=1,nnz-1
       if(cols(k).ne.cols(k+1)) then
          CSC_%colptr(n)=k+1
          n=n+1
       end if
    end do
    CSC_%colptr(CSC_%nmat+1)=nnz+1

    do k = 1, CSC_%nmat
       CSC_%loc2glb(k) = k
    enddo

    nnz = CSC_%colptr(CSC_%nmat+1)-1
    flagcor = API_YES
    call pastix_fortran_checkMatrix(CSC_%data_check, CSC_%pastix_comm, CSC_%iparm(IPARM_VERBOSE),&
         CSC_%iparm(IPARM_SYM), flagcor, CSC_%nmat, CSC_%colptr, CSC_%row, CSC_%avals,&
         CSC_%loc2glb, CSC_%iparm(IPARM_DOF_NBR))
    if (nnz.NE.(CSC_%colptr(CSC_%nmat+1)-1)) then
       deallocate(CSC_%row,CSC_%avals)
       CSC_%nnz = CSC_%colptr(CSC_%nmat+1)-1
       allocate(CSC_%row(CSC_%nnz))
       allocate(CSC_%avals(CSC_%nnz))
       call pastix_fortran_checkMatrix_End(CSC_%data_check, CSC_%iparm(IPARM_VERBOSE), &
            CSC_%row, CSC_%avals, CSC_%iparm(IPARM_DOF_NBR))
    endif
    CSC_%iparm(IPARM_MATRIX_VERIFICATION) = API_NO
  end subroutine fill_implicit_pastix


  subroutine init_implicit_pastix(CSC_)
    use all_variables, only : global_parameters, zones
    implicit none
    Type(CSC) :: CSC_
    integer*4 k,nmat
    real*8,allocatable,target :: dummy(:)
    integer :: OMP_GET_MAX_THREADS
    nmat=0
    do k=1,global_parameters%N_Zones
       nmat=nmat+(zones(k)%mesh%Nx+2)*(zones(k)%mesh%Nz+2)
    end do
    allocate(CSC_%colptr(1:nmat+1))
    allocate(CSC_%b(1:nmat))
    allocate(CSC_%loc2glb(1:nmat))
    allocate(CSC_%perm(1:nmat))
    allocate(CSC_%invp(1:nmat))

    allocate(dummy(1:1))
    CSC_%row => CSC_%colptr
    CSC_%avals => dummy

    CSC_%pastix_comm = 0  ! no mpi
    
    CSC_%nmat=nmat
    CSC_%pastix_data = 0
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

  subroutine analyze_implicit_pastix(CSC_)
    implicit none
    Type(CSC) :: CSC_
    CSC_%iparm(IPARM_START_TASK) = API_TASK_ORDERING
    CSC_%iparm(IPARM_END_TASK) = API_TASK_NUMFACT
    call pastix_fortran(CSC_%pastix_data,CSC_%pastix_comm, &
         CSC_%nmat, CSC_%colptr, CSC_%row, CSC_%avals, &
         CSC_%perm, CSC_%invp, CSC_%b, &
         CSC_%rhs, CSC_%iparm, CSC_%dparm)
  end subroutine analyze_implicit_pastix


  subroutine solve_implicit_pastix(CSC_)
    implicit none
    Type(CSC) :: CSC_
    integer*4 k
    CSC_%iparm(IPARM_START_TASK)       = API_TASK_SOLVE
    CSC_%iparm(IPARM_END_TASK)         = API_TASK_SOLVE
!!$    open(unit=101,file='RHS',status='unknown')
!!$    do k=1,CSC_%nmat
!!$       write(101,*) CSC_%b(k)
!!$    end do
!!$    close(101)
    call pastix_fortran(CSC_%pastix_data,CSC_%pastix_comm, &
         CSC_%nmat, CSC_%colptr, CSC_%row, CSC_%avals, &
         CSC_%perm, CSC_%invp, CSC_%b, &
         CSC_%rhs, CSC_%iparm, CSC_%dparm)
  end subroutine solve_implicit_pastix


end module Mpastix_solve
