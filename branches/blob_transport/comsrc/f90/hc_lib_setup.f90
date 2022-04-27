! -*-Mode:f90-*-
! Hydrocarbon_Library_Setup.f90
! Hydrocarbon Calculation Library Setup File
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! October, 1999
!
! This module sets up each of the common blocks
! and assigns data for state-specific information.
!

Module HC_Lib_Setup

  Use ComHC
  Implicit None

  ! Keep type declarations in memory.
  Save

  ! Set all integer parameter variables.
  Integer, Parameter :: Max_Reactions_Per_State = 50
  Integer, Parameter :: State_Name_Length = 10
  ! jdemod - include entire reaction description as part of the reaction table
  integer, parameter :: Max_reaction_string_len=100

  ! Declare state entry Table format for HCs and Hs (H/D/T).
  Type State_Entry
     Integer :: State_Number
     Character (State_Name_Length) :: State_Name
     Integer :: Carbon_Content
     Integer :: Hydrogen_Content
     Integer :: State_Charge
     Real    :: State_Mass

     !
     ! jdemod - add the list of reactions leading from this state to the state information
     !          if the number_reactions is zero and the state is not the end state of the 
     !          reaction chain then there is a problem with the input atomic physics data
     !          since there is a deadend channel in the data
     !        - NOTE : can't be done this way because state_entry is used for both HC and H states
     !
     !Integer :: Number_Reactions
     !Integer, Dimension (Max_Reactions_Per_State) :: Reaction

  End Type State_Entry

  ! Define the reaction table to store state numbers, number of reactions each can
  ! undergo and the reaction numbers for each transition.  Currently, it has a 
  ! maximum of 50 reactions per state.
  Type Reaction_Table_Type
     Integer :: State_Number
     Integer :: Number_Reactions
     Integer, Dimension (Max_Reactions_Per_State) :: Reaction
  End Type Reaction_Table_Type


  ! Define the state transform table to store start and end states and whether 
  ! it undergoes a proton (p) or electron (e) reaction to make the state change.
  ! Additionally, store bounds of applicability, error, and reaction energetics.
  Type State_Transform_Table_Type
     Integer :: Reaction_Number
     Integer :: Start_State
     Integer, Dimension (Highest_Carbon_Content) :: End_C_States
     Integer, Dimension (Number_H_Products)  :: End_H_States
     Character (len = 10) :: Reaction_Identifier
     Character (len = 1) :: Reaction_Type
     character (len = max_reaction_string_len) :: reaction_desc
     real :: reaction_product_mass
     real :: charged_reaction_product_mass
     Integer :: HC_E_Type
     Real :: HC_E
     Integer :: H_E_Type
     Real :: H_E
     Integer :: Sigma_TPD
     Integer :: Sigma_Polynomial_Terms
     Integer :: SigmaV_TPD
     Integer :: SigmaV_Polynomial_Terms
     Real :: Sigma_Tmin_Limit
     Real :: Sigma_ValueEMin_Limit
     Real :: Sigma_ValueEMax_Limit
     Real :: Sigma_Energy_Error_Limit
     Real :: SigmaV_Tmin_Limit
     Real :: SigmaV_ValueEMin_Limit
     Real :: SigmaV_ValueEMax_Limit
     Real :: SigmaV_Energy_Error_Limit

     !
     ! jdemod - include a pointer to the tabulated sigma-v data for this reaction - storage will be allocated
     !          for electron reactions this is indexed by the binned electron temperature
     !          for proton reactions this is indexed by the binned electron temperatures and impacting proton energies
     !
     ! WARNING NOTE: Some of the data in the JR tables requires double precision to represent numerically - the 
     !               exponents can be as small as -160 or smaller. The PGI compiler maps these values to zero
     !               when performing the type conversion - this is fine and will work - however, other compilers
     !               might map the data differently and so this should be looked at when porting the code. 
     !
     real,pointer, dimension(:) :: t_index
     real,pointer, dimension(:) :: e_index
     real,pointer, dimension(:,:) :: reaction_data

  End Type State_Transform_Table_Type


  ! Define the re-emissions table to store the most probably new state that a
  ! hydrocarbon striking the wall will be re-released as if it does not get
  ! stuck, ie. it reflects physically or is sputtered chemically.  
  ! Note support for up to 3 re-emitted hydrocarbon species upon a chemical
  ! release, same as upon breakup from an in-flight collision.
  Type Reflection_Table_Type
     Integer :: State_Number
     Integer :: Reflected_State_Number
  End Type Reflection_Table_Type
  Type Sputtering_Table_Type
     Integer :: State_Number
     Integer, Dimension (Highest_Carbon_Content) :: Sputtered_State_Number
  End Type Sputtering_Table_Type

End Module HC_Lib_Setup
