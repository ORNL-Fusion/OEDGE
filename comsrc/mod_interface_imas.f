!     -*-Fortran-*-
!
!     ==================================================================
!     ==================================================================
!
!      MODULE MOD_INTERFACE_TAGS
!      IMPLICIT none
!      PUBLIC

!...  Data types: 
!      INTEGER, PARAMETER, PUBLIC :: DTY_B = 1,     ! byte
!     .                              DTY_I = 2,     ! integer (4 byte)
!     .                              DTY_R = 3,     ! single precision real (4 byte)
!     .                              DTY_D = 4      ! double precision real (8 byte)    

!      END MODULE MOD_INTERFACE_TAGS
!
!     ==================================================================
!     ==================================================================
!
      MODULE MOD_INTERFACE_IMAS
      USE ids_schemas
      USE ids_routines
      IMPLICIT none
      PRIVATE


!      INTERFACE inPutData
!        MODULE PROCEDURE PutDataI ,PutDataR ,PutDataD,
!     .                   PutDatumI,PutDatumR,PutDatumD
!      END INTERFACE

!      INTERFACE inGetData
!        MODULE PROCEDURE GetDataI ,GetDataR ,GetDataD,
!     .                   GetDatumI,GetDatumR,GetDatumD
!      END INTERFACE


!...  Definitions:
!     ==================================================================

!
!     ------------------------------------------------------------------
!      TYPE :: type_interface
!        REAL           :: version
!        CHARACTER(512) :: file_name
!        INTEGER        :: file_pointer
!        LOGICAL        :: file_open
!        LOGICAL        :: file_stream
!      ENDTYPE type_interface
!
!     ------------------------------------------------------------------


!..  .Declarations:
!     ==================================================================

      PUBLIC :: imasTest

!      PUBLIC :: inOpenInterface, inCloseInterface, inPutData , inGetData

!      TYPE(type_interface) :: interface


!...  Routines:
!     ==================================================================
!
      CONTAINS


      SUBROUTINE imasTest
      USE mod_sol28_global
      IMPLICIT none

!     INTEGER  , INTENT(IN) :: io_select
!     CHARACTER, INTENT(IN) :: file_name*(*)

      if (opt%log.GE.2) write(0,*) 'here in imas-land'

      RETURN
 99   STOP
      END SUBROUTINE imasTest

!
!     ---------------------------------------------------------------------
!     Public: 
!     ---------------------------------------------------------------------
!
!      SUBROUTINE inOpenInterface(file_name,io_select)
!      IMPLICIT none
!
!      INTEGER  , INTENT(IN) :: io_select
!      CHARACTER, INTENT(IN) :: file_name*(*)
!
!      RETURN
! 99   STOP
!      END SUBROUTINE inOpenInterface
!
!     ---------------------------------------------------------------------
!
      END MODULE MOD_INTERFACE_IMAS


