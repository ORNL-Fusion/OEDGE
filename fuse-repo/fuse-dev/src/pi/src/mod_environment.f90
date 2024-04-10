! ==================================================================
! ==================================================================
!
!MODULE MOD_place_holder
!  IMPLICIT none
!  PUBLIC
!
! Data types: 
!  INTEGER, PARAMETER, PUBLIC :: DTY_B = 1,&    ! byte
!                                DTY_I = 2,&    ! integer (4 byte)
!                                DTY_R = 3,&    ! single precision real (4 byte)
!                                DTY_D = 4      ! double precision real (8 byte)    
!
!END MODULE MOD_place_holder
!
! ==================================================================
! ==================================================================
!
MODULE mod_environment
#include "compile_opt.inc"
!     USE mod_interface_tags
  IMPLICIT none
  PRIVATE

!  INTERFACE inPutData
!    MODULE PROCEDURE PutDataI ,PutDataR ,PutDataD,
!                     PutDatumI,PutDatumR,PutDatumD
!  END INTERFACE

!
! ------------------------------------------------------------------

  TYPE :: type_interface
    REAL           :: version
    CHARACTER(512) :: file_name
    INTEGER        :: file_pointer
    LOGICAL        :: file_open
    LOGICAL        :: file_stream
  ENDTYPE type_interface

  INTEGER :: MAXNDAT = 300

! Declarations:
! ==================================================================

  PUBLIC :: piInvestigateEnvironment

  TYPE(type_interface) :: interface
     
  INTEGER :: ndat

! Data types: 
  INTEGER, PARAMETER, PUBLIC :: ITF_READ  = 1  ! open interface for reading data

! Routines:
! ==================================================================
!
  CONTAINS
!
!   ---------------------------------------------------------------------
!   Private: 
!   ---------------------------------------------------------------------
!
  SUBROUTINE piInvestigateCudaFortran
#if CUDA_FORTRAN == 1
  USE cudafor
  IMPLICIT none

  type (cudaDeviceProp) :: prop
  integer :: nDevices=0, i, ierr

  ierr = cudaGetDeviceCount(nDevices)

  do i = 0, nDevices-1
    write(*,"('Device Number: ',i0)") i
    ierr = cudaGetDeviceProperties(prop, i)
    write(*,"(' Device Name: ',a)") trim(prop%name)
    write(*,"(' Compute Capability: ',i0,'.',i0)") &
    prop%major, prop%minor
    write(*,"(' Number of Multiprocessors: ',i0)") &
    prop%multiProcessorCount
    write(*,"(' Max Threads per Multiprocessor: ',i0)") &
    prop%maxThreadsPerMultiprocessor
    write(*,"(' Global Memory (GB): ',f9.3,/)") &
    prop%totalGlobalMem/1024.0**3
  end do

    
# endif      
  RETURN
  END SUBROUTINE piInvestigateCudaFortran
!
!   ---------------------------------------------------------------------
!
  SUBROUTINE piInvestigateMPI
  IMPLICIT none
#if MPI == 1
  include 'mpif.h'

  INTEGER :: numtasks, rank, len, ierr  
  CHARACTER(MPI_MAX_PROCESSOR_NAME) :: hostname
  
  WRITE(0,*) 'investigating MPI' 

! from https://computing.llnl.gov/tutorials/mpi/ 
      
  ! initialize MPI
  call MPI_INIT(ierr)

  ! get number of tasks
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

  ! get my rank
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

  ! this one is obvious
  call MPI_GET_PROCESSOR_NAME(hostname, len, ierr)
  print *, 'Number of tasks=',numtasks,' My rank=',rank,' Running on=',hostname

  ! do some work with message passing 

  ! done with MPI
  call MPI_FINALIZE(ierr)
# endif      
  RETURN
  END SUBROUTINE piInvestigateMPI
!
! ---------------------------------------------------------------------
! Public:
! ---------------------------------------------------------------------
!
  SUBROUTINE piInvestigateEnvironment
  IMPLICIT none

#if MPI == 1
  CALL piInvestigateMPI
#endif

#if CUDA_FORTRAN == 1
  CALL piInvestigateCudaFortran
#endif
      
  RETURN
  END SUBROUTINE piInvestigateEnvironment
!     
!      ==================================================================
!
END MODULE mod_environment


