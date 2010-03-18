! -*-Mode:f90-*-
! Initialize_OUT_Data.f90
! OUT Program Data Storage Initialization File
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! December, 2002
!
! This module initializes each of the required data structures for
! the OUT program.  It acts as the data-passing back-end between
! the DIVIMP code and the output processing OUT routines added.

Module HC_Init_OUT_Data

  ! Use statements.
  Use HC_OUT_Storage_Setup ! Access to derived data types.	

  ! Every good Fortran 90 program has...
  Implicit None

  ! Define data tables.
  !Type (HC_OUT_Data_Table_Type) :: HC_OUT_Data_Table

Contains

  Subroutine Initialize_HC_OUT_Data_Table

    ! Every good Fortran 90 program has...
    Implicit None

    HC_State_List = ""
    HC_Density = 0.0
    H_Density = 0.0
    HC_Output_List = 0
    HC_Walks = 0.0
    HC_Factor_A = 0.0
    HC_Factor_B = 0.0

    HC_Trans_Prob = 0.0
    HC_KVALS = 0.0
    HC_FP = 0.0
    HC_FT = 0.0

  End Subroutine Initialize_HC_OUT_Data_Table

End Module HC_Init_OUT_Data
