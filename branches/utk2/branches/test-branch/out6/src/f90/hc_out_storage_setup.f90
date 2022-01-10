! -*-Mode:f90-*-
! Storage_setup.f90
! Storage Setup File
!
! Adam McLean under Prof. Peter Stangeby
! Research Assistant: David Elder
! January, 2002
!
! This module sets up each of the common blocks
! and assigns data for cell-specific information.
!

Module HC_OUT_Storage_Setup

  Use ComHC ! Access to setup and input values.

  Implicit None

  ! Keep type declarations in memory.
  Save



  !Type HC_OUT_Data_Table_Type

  ! Contains all data read into OUT from the DIV binary output file.
  Character (Len=10), Dimension (Number_HC_Species) :: HC_State_List ! List of character names of each hydrocarbon species used in OUT processing.
  Real, Dimension (maxnks,maxnrs,Number_HC_Species) :: HC_Density ! Hydrocarbon density array. Records cell position and species of each fragment after each timestep.
  Real, Dimension (maxnks,maxnrs,Num_H_States) :: H_Density ! Hydrogen density array. Records cell position and species of each fragment at time of discard during breakup.
  Integer, Dimension (Number_HC_Species+1) :: HC_Output_List ! Determines which hydrocarbons to include in stored binary file for OUT processing: 0=not included, 1=included.
  Real, Dimension (Max_Number_Walks,2) :: HC_Walks ! Follows position of particles as they travel for each timestep up to maxnws steps - typically 10,000.
  ! jdemod - change dimensions to appropriate values
  Real, Dimension (Number_HC_Species) :: HC_Factor_A ! Stores FACTA.
  Real, Dimension (Number_HC_Species) :: HC_Factor_B ! Stores FACTB.
  !Real, Dimension (-1:Number_HC_Species) :: HC_Factor_A ! Stores FACTA.
  !Real, Dimension (-1:Number_HC_Species) :: HC_Factor_B ! Stores FACTB.
  Real, Dimension (maxnks,maxnrs,-1:0) :: HC_TIZS_CH ! CH ionization locations. (:,:,-1) for primary neutrals and (:,:,0) for those from primaries + reflections/self-sputters.
  !		Real, Dimension (maxnks,maxnrs,-1:0) :: HC_TIZS_C2 ! C2 ionization locations. (:,:,-1) for primary neutrals and (:,:,0) for those from primaries + reflections/self-sputters.

  ! Storage used solely in OUT.
  Integer :: HC_Particles_Launched ! Number of hydrocarbons launched.
  Real, Dimension (maxnks,maxnrs) :: HC_Trans_Prob ! Stores total transition probability to move out of one state given reactions possible and cell properties.
  Real, Dimension (maxnks,Number_HC_Species) :: HC_KVALS
  Real, Dimension (maxnks,maxnrs) :: X_CH ! CH/CD excitation rate.
  Real, Dimension (maxnks,maxnrs) :: X_C2 ! C2 excitation rate.
  Real :: HC_FP ! Density normalization factor.
  Real :: HC_FT ! Density normalization factor.
  Real :: HC_ABSFAC ! Normalization for wall launches.
  !		Real :: FYTOT

  !End Type HC_Out_Data_Table_Type

End Module HC_OUT_Storage_Setup
